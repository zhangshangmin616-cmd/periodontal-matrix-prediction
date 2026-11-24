"""Locate and retrieve the NCBI-issued uniform counts for the designated culture libraries."""
import csv
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/"stage5/zhao"))
import download_ncbi_counts as ncbi

RUNS={"SRR6291210","SRR6291212"}


def main():
    a,b=0,ncbi.byte_slice(0,300)[1]
    while b-a>65536:
        mid=(a+b)//2
        data,_=ncbi.byte_slice(mid,mid+1000)
        if data.splitlines()[1].split(b"\t")[0].decode()<min(RUNS):
            a=mid
        else:
            b=mid
    data,_=ncbi.byte_slice(max(0,a-65536),b+131072)
    rows=[line.decode().split("\t") for line in data.splitlines()[1:-1]]
    assert rows[0][0]<min(RUNS) and rows[-1][0]>max(RUNS)
    found=[dict(zip(ncbi.FIELDS,r,strict=True)) for r in rows if r[0] in RUNS and r[1]=="PRJNA418326"]
    with (HERE/"count_manifest.tsv").open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=ncbi.FIELDS,delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(found)
    print(len(found),"NCBI run-level counts located",flush=True)
    assert {r["run"] for r in found}==RUNS
    ncbi.ROOT=HERE
    (HERE/"counts").mkdir(exist_ok=True)
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(ncbi.fetch,found))


if __name__=="__main__":
    main()
