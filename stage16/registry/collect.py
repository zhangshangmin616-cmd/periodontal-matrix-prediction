"""Retrieve public discovery metadata. Run with --out to preserve the archived snapshot."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from html import unescape
import io
import json
from pathlib import Path
import re
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET


DATA = "https://api.datacite.org/dois"
ENTREZ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
BIO = "https://www.ebi.ac.uk/biostudies/api/v1/search"
ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
SRA_TERMS = [
    'txid9606[Organism] AND (orthodontic OR "tooth movement")',
    'txid9606[Organism] AND ("periodontal ligament" OR "periodontal membrane") AND ("pressure side" OR "tension side" OR "compression side" OR "pressure-side" OR "tension-side")',
    'txid9606[Organism] AND ("periodontal ligament" OR "periodontal membrane") AND ("7 days" OR "day 7" OR "7d" OR "one week")',
]
DC_TERMS = [
    'types.resourceTypeGeneral:Dataset AND titles.title:("orthodontic tooth movement" OR "orthodontic force" OR "tooth movement")',
    'types.resourceTypeGeneral:Dataset AND ("periodontal ligament" AND "orthodontic tooth movement")',
    'types.resourceTypeGeneral:Dataset AND ("periodontal ligament" AND (tension OR "pressure side" OR "compression side"))',
    '("pressure side" OR "compression side") AND "tension side" AND "periodontal ligament"',
    'relatedIdentifiers.relatedIdentifier:("10.1093/ejo/cjag053" OR "10.1007/s00784-021-04178-y")',
]
BS_TERMS = ['"pressure side" AND "tension side" AND "periodontal ligament"',
            '"compression side" AND "tension side" AND "periodontal ligament"']
FIELDS = "study_accession,experiment_accession,sample_accession,sample_alias,sample_title,experiment_title,run_accession,library_strategy"


def request(endpoint, params):
    address = endpoint + "?" + urlencode(params)
    for attempt in range(3):
        try:
            with urlopen(Request(address,headers={"User-Agent":"PDL-clinical-metadata-audit/1.0"}),timeout=40) as f:
                return address, f.read().decode("utf8")
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def tidy(text):
    return " ".join(unescape(re.sub("<[^>]*>", " ", text or "")).split())


def sra():
    log, ids = [], set()
    for term in SRA_TERMS:
        url, raw = request(ENTREZ+"esearch.fcgi",dict(db="sra",term=term,retmax=10000,retmode="json"))
        query = json.loads(raw)["esearchresult"]
        assert len(query["idlist"])==int(query["count"])
        ids.update(query["idlist"])
        log.append(dict(query=term,url=url,count=int(query["count"]),uids=query["idlist"],
                        translation=query["querytranslation"],warnings=query.get("warninglist",{})))
        time.sleep(.35)
    ids=sorted(ids,key=int)
    summaries=[]
    for start in range(0,len(ids),40):
        batch=ids[start:start+40]
        _, raw=request(ENTREZ+"esummary.fcgi",dict(db="sra",id=",".join(batch),retmode="json"))
        response=json.loads(raw)["result"]
        assert set(response["uids"])==set(batch)
        for uid in batch:
            record=response[uid]
            meta=ET.fromstring("<root>"+record["expxml"]+"</root>")
            text=lambda path: meta.findtext(path) or ""
            study,experiment=meta.find("Study"),meta.find("Experiment")
            summaries.append(dict(uid=uid,project=text("Bioproject"),study=study.get("acc"),
                study_title=study.get("name"),experiment=experiment.get("acc"),title=text("Summary/Title"),
                biosample=text("Biosample"),strategy=text(".//LIBRARY_STRATEGY"),source=text(".//LIBRARY_SOURCE"),
                runs=re.findall('acc="([^"]+)"',record["runs"])))
        time.sleep(.35)
    return "sra.json",dict(queries=log,records=summaries)


def datacite():
    log, records=[],{}
    for term in DC_TERMS:
        page,found,urls=1,[],[]
        while True:
            params={"query":term,"page[size]":100,"page[number]":page}
            url,raw=request(DATA,params)
            data=json.loads(raw)
            urls.append(url)
            for record in data["data"]:
                a=record["attributes"]
                doi=record["id"]
                found.append(doi)
                records[doi]=dict(title=" | ".join(tidy(t["title"]) for t in a["titles"]),
                    publisher=a["publisher"],url=a["url"],type=a["types"]["resourceTypeGeneral"],
                    description=" ".join(tidy(d.get("description","")) for d in a["descriptions"]),
                    related=[{k:v for k,v in r.items() if k in ["relatedIdentifier","relationType"]}
                             for r in a.get("relatedIdentifiers",[])],registered=a["registered"])
            if len(found)>=data["meta"]["total"]:
                break
            page+=1
        assert len(set(found))==data["meta"]["total"]
        log.append(dict(query=term,urls=urls,count=len(found),doi=found))
    return "datacite.json",dict(queries=log,records=records)


def biostudies():
    log=[]
    for term in BS_TERMS:
        url,raw=request(BIO,dict(query=term,pageSize=200))
        data=json.loads(raw)
        hits=[dict(accession=v["accession"],title=tidy(v["title"])) for v in data["hits"]]
        assert len(hits)==int(data["totalHits"])
        log.append(dict(query=term,url=url,count=len(hits),hits=hits))
    return "biostudies.json",dict(queries=log)


def mirror():
    series="GSE327337"
    params=dict(acc=series,targ="self",form="text",view="full")
    geo,raw=request(GEO,params)
    samples=re.findall(r"^!Series_sample_id = (GSM\d+)",raw,re.M)
    assert len(samples)==5 and "PRJNA1450544" in raw
    url,raw=request(ENA,dict(accession="PRJNA1450544",result="read_run",fields=FIELDS,format="tsv"))
    runs=list(csv.DictReader(io.StringIO(raw),delimiter="\t"))
    assert set(samples)=={r["sample_alias"] for r in runs}
    for r in runs:
        _,raw=request(GEO,dict(acc=r["sample_alias"],targ="self",form="text",view="full"))
        r["geo_title"]=re.search(r"^!Sample_title = (.+)",raw,re.M)[1].strip()
        r["sample_relations"]=re.findall(r"^!Sample_relation = ([^\r\n]+)",raw,re.M)
        assert r["sample_title"]==r["geo_title"] and r["sample_accession"] in raw and r["experiment_accession"] in raw
    return "source_crosswalk.json",dict(geo=geo,ena=url,runs=runs)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--resources",nargs="+",choices=["sra","datacite","biostudies","mirror"],
                   default=["sra","datacite","biostudies","mirror"])
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    now=datetime.now(timezone.utc).isoformat(timespec="seconds")
    with ThreadPoolExecutor(4) as pool:
        tasks=[pool.submit(f) for f in (sra,datacite,biostudies,mirror) if f.__name__ in args.resources]
        for task in tasks:
            name,data=task.result()
            (args.out/name).write_text(json.dumps(dict(retrieval_utc=now,**data),ensure_ascii=False,indent=2)+"\n",encoding="utf8")
            print(name,flush=True)


if __name__=="__main__":
    main()
