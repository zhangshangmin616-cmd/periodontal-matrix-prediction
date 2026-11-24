"""Locate eight public SRA gene-count objects and verify their source checksums."""
from concurrent.futures import ThreadPoolExecutor
import csv
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]/'stage5/zhao'))
import download_ncbi_counts as ncbi


def main():
    ncbi.ROOT=HERE
    meta=json.loads((HERE/'ena_runs.json').read_text(encoding='utf8'))
    runs=sorted(r['run_accession'] for r in meta)
    assert len(runs)==8
    _,end=ncbi.byte_slice(0,300)
    start=0
    while end-start>65536:
        mid=(start+end)//2
        part,_=ncbi.byte_slice(mid,mid+1000)
        sample=part.splitlines()[1].split(b'\t')[0].decode()
        if sample<runs[0]:
            start=mid
        else:
            end=mid
    obj,_=ncbi.byte_slice(max(start-32768,0),end+65536)
    near=[row.decode().split('\t') for row in obj.splitlines()[1:-1]]
    matches=[dict(zip(ncbi.FIELDS,r)) for r in near if r[0] in runs and r[1]=='PRJNA665587']
    assert sorted(r['run'] for r in matches)==runs, (len(matches),runs)
    manifest=HERE/'count_manifest.tsv'
    if manifest.exists():
        with manifest.open(encoding='utf8') as f:
            assert matches==list(csv.DictReader(f,delimiter='\t'))
    else:
        with manifest.open('w',encoding='utf8',newline='') as f:
            out=csv.DictWriter(f,fieldnames=ncbi.FIELDS,delimiter='\t',lineterminator='\n')
            out.writeheader()
            out.writerows(matches)
    (HERE/'counts').mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(ncbi.fetch,matches))
    print('Verified',len(matches),'count files against original MD5.')


if __name__=='__main__':
    main()
