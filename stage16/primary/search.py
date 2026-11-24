"""Search frozen-panel names in PubMed and retain the exact requests and returned IDs."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import time
from urllib.parse import urlencode

HERE = Path(__file__).resolve().parent
TERM = ('(orthodontic*[tiab] OR "tooth movement"[tiab] OR "maxillary expansion"[tiab]) AND '
        '("periodontal ligament"[tiab] OR "periodontal tissue"[tiab] OR "periodontal tissues"[tiab]) AND '
        '(premolar*[tiab] OR patient*[tiab] OR volunteer*[tiab] OR "split-mouth"[tiab] OR '
        '"compression side"[tiab] OR "pressure side"[tiab] OR "tension side"[tiab] OR "tension sides"[tiab])')


def read(path):
    with path.open(encoding="utf8", newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def save(path,records):
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def run(helper,payload):
    for attempt in range(3):
        p=subprocess.run([sys.executable,str(helper)],input=json.dumps(payload),capture_output=True,text=True,encoding="utf8")
        res=json.loads(p.stdout)
        if res.get("ok"):
            time.sleep(.37)
            return res
        time.sleep(1+attempt)
    raise ValueError(res)


def collect(helper):
    panel=read(HERE/"query_terms.tsv")
    queries=[]
    for i in range(0,len(panel),6):
        batch=panel[i:i+6]
        phrases=sorted(set(x for g in batch for x in g["terms"].split("|")))
        names=" OR ".join('"'+x+'"[tiab]' for x in phrases)
        search=f"({names}) AND {TERM}"
        res=run(helper,{"endpoint":"esearch","params":{"db":"pubmed","term":search,"retmode":"json","retmax":500},"max_items":500})
        ids=res.get("records",[])
        assert len(ids)<500 and not res.get("truncated")
        queries.append(dict(batch=i//6+1,genes="|".join(g["gene"] for g in batch),matches=len(ids),
                            pmids="|".join(ids),query=search,
                            url="https://pubmed.ncbi.nlm.nih.gov/?"+urlencode({"term":search})))
    save(HERE/"pubmed_queries.tsv",queries)
    print("Per-query counts",[r["matches"] for r in queries])
    return sorted({pm for r in queries for pm in r["pmids"].split("|") if pm},key=int)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("entrez_helper",type=Path)
    args=cli.parse_args()
    ids=collect(args.entrez_helper)
    print("Unique PMIDs",len(ids))
