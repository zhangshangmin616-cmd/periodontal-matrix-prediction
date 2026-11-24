"""Fetch NCBI-authored counts, validate their original MD5, and store compressed copies.

No sample expression measurements are parsed during download. The original
run-index lines originate in the NCBI daily S3 content.tsv manifest.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import os
from pathlib import Path
import time
import urllib.request

ROOT=Path(__file__).resolve().parent
BUCKET="https://sra-rnaseq-analysis.s3.amazonaws.com/"
FIELDS=["run","bioproject","biosample","analysis","organism","path","md5"]
PROJECT="PRJNA1136245"
FIRST,LAST="SRR29839648","SRR29839710"
INDEX=ROOT/"ncbi_run_manifest.tsv"


def byte_slice(start,end):
    path=BUCKET+f"content.tsv?bytes={start}-{end}"
    req=urllib.request.Request(path,headers={"Range":f"bytes={start}-{end}"})
    with urllib.request.urlopen(req,timeout=40) as response:
        if response.status!=206:
            raise ValueError("Server did not honor byte range")
        length=int(response.headers["Content-Range"].rsplit("/",1)[1])
        return response.read(),length


def discover():
    _,end=byte_slice(0,300)
    start=0
    # content.tsv is lexicographically ordered on run accession. Narrow around
    # the published series, then retrieve every row in its complete interval.
    while end-start>65536:
        mid=(start+end)//2
        data,_=byte_slice(mid,mid+800)
        accession=data.splitlines()[1].split(b"\t")[0].decode()
        if accession < FIRST:
            start=mid
        else:
            end=mid
    data,_=byte_slice(max(0,start-32768),end+65536)
    rows=[line.decode().split("\t") for line in data.splitlines()[1:-1]]
    if not (rows[0][0] < FIRST and rows[-1][0] > LAST):
        raise ValueError("Read interval did not span all 63 runs")
    matches=[row for row in rows if row[1]==PROJECT]
    expected={f"SRR{i}" for i in range(int(FIRST[3:]),int(LAST[3:])+1)}
    assert len(matches)==63 and {r[0] for r in matches}==expected
    if INDEX.exists():
        with INDEX.open(encoding="utf-8") as stream:
            saved=list(csv.reader(stream,delimiter="\t"))
        if saved[1:]!=matches:
            raise ValueError("Remote manifest versions differ from frozen study copy")
    else:
        with INDEX.open("w",encoding="utf-8",newline="") as stream:
            writer=csv.writer(stream,delimiter="\t",lineterminator="\n")
            writer.writerow(FIELDS)
            writer.writerows(matches)


def fetch(rec):
    output=ROOT/"counts"/(rec["run"]+".tsv.gz")
    if output.exists():
        with gzip.open(output,"rb") as stream:
            if hashlib.file_digest(stream,"md5").hexdigest()==rec["md5"]:
                return
    tmp=output.with_suffix(".part")
    for attempt in range(4):
        try:
            check=hashlib.md5()
            with urllib.request.urlopen(BUCKET+rec["path"],timeout=50) as response, tmp.open("wb") as file:
                with gzip.GzipFile(filename="",mode="wb",fileobj=file,mtime=0,compresslevel=6) as compressed:
                    while chunk:=response.read(65536):
                        check.update(chunk)
                        compressed.write(chunk)
            if check.hexdigest()!=rec["md5"]:
                raise ValueError(f"Checksum mismatch {rec['run']}")
            os.replace(tmp,output)
            return
        except (OSError,ValueError):
            if attempt==3:
                raise
            time.sleep(1+attempt)


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("--discover",action="store_true",help="Verify manifest via HTTP Range on daily sorted run index")
    cmd=args.parse_args()
    if cmd.discover or not INDEX.exists():
        discover()
    (ROOT/"counts").mkdir(exist_ok=True)
    with INDEX.open(encoding="utf-8") as stream:
        records=list(csv.DictReader(stream,delimiter="\t"))
    with ThreadPoolExecutor(max_workers=10) as pool:
        list(pool.map(fetch,records))
    print(f"{len(records)} gene-count files stored with original MD5 verified.")
