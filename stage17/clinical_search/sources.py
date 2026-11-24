"""Re-fetch original public files or check a few read records per technical library."""
import argparse
import csv
import hashlib
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request,urlopen

HERE = Path(__file__).resolve().parent
FILES = {
    "series.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189288/matrix/GSE189288_series_matrix.txt.gz",
    "counts.gz": "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE189288&format=file&file=GSE189288_raw_counts_GRCh38.p13_NCBI.tsv.gz",
}
for project in ["PRJNA782351","PRJNA698902"]:
    FILES[project+".tsv"]=("https://www.ebi.ac.uk/ena/portal/api/filereport?accession="+project+
        "&result=read_run&fields=study_accession,sample_accession,run_accession,experiment_accession,"
        "sample_alias,experiment_title,library_name,fastq_md5,read_count,base_count,fastq_ftp&format=tsv")


def prefix(r):
    request=Request("https://"+r["fastq_ftp"].split(";")[0],headers={"Range":"bytes=0-32767","User-Agent":"PDL-provenance/1.0"})
    with urlopen(request,timeout=20) as stream:
        first=stream.read(32768)
    fields=zlib.decompressobj(31).decompress(first).splitlines()
    assert fields[0].startswith(b"@") and fields[2].startswith(b"+")
    digest=hashlib.sha256(b"\n".join(fields[1:80:4])).hexdigest()
    return (r["study_accession"],r["sample_alias"],r["run_accession"],20,digest)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", action="store_true",help="Obtain fingerprints of the first 20 forward reads per library.")
    parser.add_argument("--outdir",type=Path,default=HERE)
    args=parser.parse_args()
    args.outdir.mkdir(parents=True,exist_ok=True)
    if args.reads:
        runs=[]
        for project in ["PRJNA782351","PRJNA698902"]:
            with (HERE/(project+".tsv")).open(encoding="utf8") as f:
                runs.extend(csv.DictReader(f,delimiter="\t"))
        with ThreadPoolExecutor(6) as pool:
            rows=list(pool.map(prefix,runs))
        assert len({r[-1] for r in rows})==len(rows)==51
        with (args.outdir/"fastq_prefixes.tsv").open("w",encoding="utf8",newline="") as f:
            w=csv.writer(f,delimiter="\t",lineterminator="\n")
            w.writerow(["project","sample","run","first_forward_reads","sha256_read_sequences"])
            w.writerows(sorted(rows))
    else:
        for name,uri in FILES.items():
            with urlopen(uri,timeout=45) as stream:
                data=stream.read()
            (args.outdir/name).write_bytes(data)
            print(name,len(data),hashlib.sha256(data).hexdigest())


if __name__=="__main__":
    main()
