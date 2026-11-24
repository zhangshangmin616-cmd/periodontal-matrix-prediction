"""Find public NCBI gene counts for the Wuhan control-vector libraries."""
import csv
import io
import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "stage5/zhao"))
import download_ncbi_counts as ncbi

PROJECT = "PRJNA1377575"
SOURCE = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + PROJECT +
          "&result=read_run&fields=run_accession,sample_alias,experiment_title&format=tsv")


def write_table(path, rows):
    with path.open("w", encoding="utf8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    with urlopen(SOURCE, timeout=40) as res:
        raw = list(csv.DictReader(io.StringIO(res.read().decode()), delimiter="\t"))
    selected = []
    for row in raw:
        match = re.search(r"sh-Control, tension (0|24|48) h, rep ([123])", row["experiment_title"])
        if match:
            selected.append(dict(run=row["run_accession"], gsm=row["sample_alias"], hours=match[1],
                                 replicate=match[2], experiment_title=row["experiment_title"]))
    selected.sort(key=lambda row: row["run"])
    assert len({r["gsm"] for r in selected}) == 9 and len({r["hours"]+":"+r["replicate"] for r in selected}) == 9
    write_table(HERE / "libraries.tsv", selected)
    runs = {r["run"] for r in selected}
    low, high = 0, ncbi.byte_slice(0,300)[1]
    while high - low > 65536:
        mid = (low + high) // 2
        fragment, _ = ncbi.byte_slice(mid, mid+1100)
        last = fragment.splitlines()[1].split(b"\t")[0].decode()
        if last < min(runs):
            low = mid
        else:
            high = mid
    window, _ = ncbi.byte_slice(max(low-32768,0), high+131072)
    records = [s.decode().split("\t") for s in window.splitlines()[1:-1]]
    assert records[0][0] < min(runs) and records[-1][0] > max(runs)
    found = [dict(zip(ncbi.FIELDS, r, strict=True)) for r in records if r[0] in runs and r[1] == PROJECT]
    (HERE/"availability.json").write_text(json.dumps(dict(ENA=SOURCE, runs_expected=len(runs),
        found=len(found), absent=sorted(runs-{r["run"] for r in found})), indent=2), encoding="utf8")
    print("Expected",len(runs),"NCBI count tables. Found",len(found),flush=True)
    if len(found) != len(runs):
        return
    write_table(HERE/"count_manifest.tsv", found)
    ncbi.ROOT = HERE
    (HERE/"counts").mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=5) as ex:
        list(ex.map(ncbi.fetch, found))
    print("Source MD5s matched all downloaded count files.", flush=True)


if __name__ == "__main__":
    main()
