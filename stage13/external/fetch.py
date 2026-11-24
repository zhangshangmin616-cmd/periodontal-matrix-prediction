"""Download primary-source files to the local raw directory and verify digests."""
import concurrent.futures
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
RAW=HERE/"raw"


def fetch(item):
    filename, url, expected = item
    path=RAW/filename
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==expected:
        return
    with urlopen(Request(url, headers={"User-Agent":"Mozilla/5.0"}),timeout=60) as response:
        body=response.read()
    if hashlib.sha256(body).hexdigest()!=expected:
        raise ValueError("Primary download no longer matches recorded bytes: "+url)
    path.write_bytes(body)


def main():
    RAW.mkdir(exist_ok=True)
    hk=json.loads((HERE/"source_manifest.json").read_text(encoding="utf8"))
    osaka=json.loads((HERE/"chip_qc.json").read_text(encoding="utf8"))
    digests={**hk,**{k[:-7]:v for k,v in osaka.items() if k.endswith("_sha256")}}
    jobs=[]
    for name, sha in digests.items():
        if name.startswith("GSE173891_"):
            url="https://ftp.ncbi.nih.gov/geo/series/GSE173nnn/GSE173891/suppl/"+name
        elif name.endswith(".txt.gz"):
            gsm=name.split(".")[0]
            url=f"https://ftp.ncbi.nih.gov/geo/samples/GSM474nnn/{gsm}/suppl/{name}"
        else:
            url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={name}&targ=self&form=text&view=full"
        jobs.append((name,url,sha))
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fetch,jobs))


if __name__ == "__main__":
    main()
