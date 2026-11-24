"""Extract matched external culture measurements using stable source identifiers."""
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
WORK = HERE.parents[1]


def identifier_map(rows, field):
    by_id = defaultdict(set)
    for row in rows:
        if row[field]:
            by_id[row[field]].add(row["symbol"])
    return {key: next(iter(names)) for key, names in by_id.items() if len(names) == 1}


def main():
    with gzip.open(WORK / "stage5/paired/hgnc_ids.tsv.gz", "rt", encoding="utf8") as stream:
        official = list(csv.DictReader(stream, delimiter="\t"))
    ens = identifier_map(official, "ensembl_gene_id")
    with (WORK / "stage5/reference/features.tsv").open(encoding="utf8") as stream:
        features = list(csv.DictReader(stream, delimiter="\t"))
    genes = [r["gene"] for r in features]
    is_target = {r["gene"] for r in features if r["ecm_component"] == "other"}
    counts = {}
    duplicate = set()
    columns = [f"{letter}{i}" for letter in ("A", "B") for i in range(1, 6)]
    totals = np.zeros(10, dtype=np.int64)
    inputfile = RAW / "GSE173891_raw_gene_count.txt.gz"
    with gzip.open(inputfile, "rt", encoding="utf8") as stream:
        reader = csv.reader(stream, delimiter="\t")
        assert next(reader) == columns
        for row in reader:
            values = np.array(row[1:], dtype=np.int64)
            assert len(values) == 10 and (values >= 0).all()
            totals += values
            name = ens.get(row[0].split(".")[0])
            if name in counts:
                duplicate.add(name)
            counts[name] = values
    retained = [gene for gene in genes if gene in counts and gene not in duplicate]
    sample_rows = []
    for i, col in enumerate(columns):
        gsm = f"GSM{5282420+i}"
        record = (RAW / gsm).read_text(encoding="utf8")
        group = "control" if i < 5 else "tension"
        number = i % 5 + 1
        assert f"!Sample_title = PDL_{group}_{number}" in record
        assert f"!Sample_description = {col}" in record
        sample_rows.append([gsm, col, group, str(number), str(totals[i])])
    with (HERE / "sample_map.tsv").open("w", encoding="utf8", newline="") as stream:
        w = csv.writer(stream, delimiter="\t", lineterminator="\n")
        w.writerow(["GSM", "column", "arm", "donor", "all_feature_total"])
        w.writerows(sample_rows)
    target = np.isin(retained, list(is_target))
    np.savez_compressed(HERE / "paired_counts.npz", genes=np.array(retained), counts=np.stack([counts[g] for g in retained]),
                        total=totals, target=target)
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [inputfile] + [
        RAW / row[0] for row in sample_rows]}
    (HERE / "source_manifest.json").write_text(json.dumps(hashes, indent=2) + "\n", encoding="utf8")
    print(len(retained), int(target.sum()), "retained background and target identifiers")
    print("Unmapped targets", sorted(is_target.difference(retained)))


if __name__ == "__main__":
    main()
