"""Download the bibliographic neighbourhood of specified published primary studies."""
import csv
import json
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode, quote
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
BASE="https://www.ebi.ac.uk/europepmc/webservices/rest/"
SEEDS={"Garlet2007":"10.1111/j.1600-0722.2007.00469.x","Garlet2008":"10.1016/j.cyto.2008.03.003",
       "Menezes2008":"10.1016/j.joen.2008.05.002","Anastasi2008":"PMID:18425345",
       "Militi2019":"10.1016/j.heliyon.2019.e02572","Chen2024":"10.1038/s41368-023-00268-7",
       "Zhao2026":"10.1177/00220345261491310"}


def fetch(path,params):
    u=BASE+path+"?"+urlencode(params)
    for retry in range(3):
        try:
            with urlopen(Request(u,headers={"User-Agent":"LiteratureSearch/1.0"}),timeout=25) as f:
                return json.load(f)
        except Exception:
            if retry==2:
                raise
            time.sleep(1.5)


def around(seed):
    name,doi=seed
    if doi.startswith("PMID:"):
        source,id_="MED",doi.split(":")[1]
    else:
        hits=fetch("search",{"query":'DOI:"'+doi+'"',"format":"json"})["resultList"]["result"]
        assert len(hits)==1,(doi,hits)
        source,id_=hits[0]["source"],hits[0]["id"]
    out=[]
    for method,folder,singular in [("citations","citationList","citation"),("references","referenceList","reference")]:
        page=1
        while True:
            response=fetch(f"{source}/{id_}/{method}",{"format":"json","page":page,"pageSize":1000})
            block=response.get(folder,{}).get(singular,[])
            for item in block:
                title=item.get("title") or item.get("citationOrder", "")
                out.append({"seed":name,"direction":method,"seed_id":id_,"id":item.get("id",""),"source":item.get("source",""),
                            "year":item.get("pubYear","") or item.get("year",""),"title":str(title).replace("\r"," ").replace("\n"," ")})
            if page*1000>=response.get("hitCount",0) or not block:
                break
            page+=1
    return out


def openalex(seed):
    name,accession=seed
    ep="https://api.openalex.org/works"
    def get(url,params):
        with urlopen(Request(url+"?"+urlencode(params),headers={"User-Agent":"LiteratureSearch/1.0"}),timeout=35) as f:
            return json.load(f)
    alias=accession.lower() if accession.startswith("PMID:") else "https://doi.org/"+accession
    work=get(ep+"/"+alias,{"select":"id,referenced_works"})
    wid=work["id"].rsplit("/",1)[-1]
    out=[]
    select="id,doi,display_name,publication_year"
    def add_set(q,mode):
        page=1
        while True:
            result=get(ep,{"filter":q,"per_page":100,"page":page,"select":select})
            for item in result["results"]:
                out.append(dict(seed=name,direction=mode,seed_id=wid,id=item["id"].rsplit("/",1)[-1],
                                source="OpenAlex",year=item["publication_year"],title=item["display_name"],
                                doi=item.get("doi") or ""))
            if page*100>=result["meta"]["count"]: break
            page+=1
    add_set("cites:"+wid,"citations")
    refs=[s.rsplit("/",1)[-1] for s in work["referenced_works"]]
    for i in range(0,len(refs),45):
        add_set("openalex:"+"|".join(refs[i:i+45]),"references")
    return out


if __name__=="__main__":
    import sys
    epmc="--epmc" in sys.argv
    with ThreadPoolExecutor(max_workers=3) as ex:
        data=[r for res in ex.map(around if epmc else openalex,SEEDS.items()) for r in res]
    with (HERE/("epmc_links.tsv" if epmc else "links.tsv")).open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(data)
    print(len(data),"links")
