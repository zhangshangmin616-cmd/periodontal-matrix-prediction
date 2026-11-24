"""Retrieve repository metadata used in the day-7 clinical periodontal search."""
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from itertools import product
from pathlib import Path
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
BASE = "https://ddbj.nig.ac.jp/search/api/entries/"
WORDS = ("orthodontic", '"periodontal ligament"', '"tooth movement"')
TYPES = ("gea", "jga-study", "sra-study", "bioproject")


def fetch(url):
    for attempt in range(3):
        try:
            with urlopen(Request(url, headers={"Accept": "application/json",
                                               "User-Agent":"PDL-metadata-audit/1.0"}), timeout=30) as fp:
                return json.load(fp)
        except OSError:
            if attempt == 2:
                raise
            sleep(1+attempt)


def query(spec):
    kind, word = spec
    params = {"keywords":word, "includeDbXrefs":"false", "includeFacets":"false",
              "perPage":"100", "fields":"identifier,title,description,organism"}
    items = {}
    page = 1
    while True:
        url = BASE+kind+"/?"+urlencode(dict(params,page=page))
        if page == 1:
            first = url
        response = fetch(url)
        data = response["pagination"]
        if page == 1:
            total = data["total"]
        for r in response["items"]:
            items[r["identifier"]] = {k:r.get(k) for k in ("identifier","title","description","organism")}
        if not data["hasNext"]:
            break
        page += 1
        if page > 100:
            raise RuntimeError("DDBJ maximum pagination reached")
    assert len(items)==total,(kind,word,total,len(items))
    return {"type":kind,"term":word,"first_page":first,"pages":page,"total":total,"items":list(items.values())}


def write_tsv(file, rows):
    with file.open("w",encoding="utf-8",newline="") as f:
        out=csv.DictWriter(f,fieldnames=rows[0],delimiter="\t",lineterminator="\n")
        out.writeheader()
        out.writerows(rows)


def main():
    with ThreadPoolExecutor(3) as pool:
        ans=list(pool.map(query,product(TYPES,WORDS)))
    write_tsv(ROOT/"ddbj_queries.tsv",[{k:v for k,v in a.items() if k!="items"} for a in ans])
    controls = []
    for kind in ("gea","jga-study"):
        url = BASE + kind + "/?" + urlencode({"perPage":"1","keywords":"cancer"})
        data = fetch(url)
        controls.append({"type":kind,"term":"cancer","url":url,
                         "total":data["pagination"]["total"],
                         "first_accession":data["items"][0]["identifier"]})
    write_tsv(ROOT/"ddbj_controls.tsv",controls)
    allrows={}
    for a in ans:
        for r in a["items"]:
            k=(a["type"],r["identifier"])
            if k not in allrows:
                allrows[k]={"type":a["type"],"identifier":r["identifier"],"matched_queries":[],
                    "title":r["title"],
                    "organism":json.dumps(r["organism"],ensure_ascii=False),
                    "url":f"https://ddbj.nig.ac.jp/search/entry/{k[0]}/{k[1]}"}
            allrows[k]["matched_queries"].append(a["term"])
    for r in allrows.values():
        r["matched_queries"]=" | ".join(r["matched_queries"])
        r["title"]=" ".join((r["title"] or "").split())
    write_tsv(ROOT/"ddbj_records.tsv",[allrows[k] for k in sorted(allrows)])
    print([(a["type"],a["term"],a["total"]) for a in ans])
    print(len(allrows),"unique repository records")


if __name__=="__main__":
    main()
