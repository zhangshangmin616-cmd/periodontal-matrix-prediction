"""Check the dated public resource indexes, deduplication and source crosswalk offline."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load(name):
    return json.loads((HERE/name).read_text(encoding="utf8"))


def main():
    sra = load("sra.json")
    records = {r["uid"]:r for r in sra["records"]}
    assert set(records)==set.union(*(set(q["uids"]) for q in sra["queries"]))
    projects = {r["project"] for r in records.values()}
    with (HERE/"projects.tsv").open(encoding="utf8",newline="") as f:
        dispositions=list(csv.DictReader(f,delimiter="\t"))
    assert {r["project"] for r in dispositions}==projects
    for q in sra["queries"]:
        assert len(set(q["uids"]))==q["count"]
    side = {records[x]["project"] for x in sra["queries"][1]["uids"]}
    day = {records[x]["project"] for x in sra["queries"][2]["uids"]}
    assert side=={"PRJNA1450544"} and day==side|{"PRJNA390189"}
    dc=load("datacite.json")
    ids=set.union(*(set(q["doi"]) for q in dc["queries"]))
    assert ids==set(dc["records"])
    titles={dc["records"][i]["title"] for i in ids}
    reports=load("biostudies.json")
    for q in reports["queries"]:
        assert q["count"]==len(q["hits"]) and all(r["accession"].startswith("S-EPMC") for r in q["hits"])
    archive = load("source_crosswalk.json")
    from_sra={r["experiment"]:r for r in records.values() if r["project"]=="PRJNA1450544"}
    assert len(archive["runs"])==len(from_sra)==5
    for r in archive["runs"]:
        matched=from_sra[r["experiment_accession"]]
        assert r["sample_alias"] in matched["title"] and r["sample_accession"]==matched["biosample"]
        assert r["sample_title"]==r["geo_title"] and matched["runs"]==[r["run_accession"]]
    assert {(r["sample_alias"],r["run_accession"]) for r in archive["runs"] if "day 7" in r["geo_title"]}=={
        ("GSM9653910","SRR37998087"),("GSM9653912","SRR37998085")}
    print(f"Checked {len(records)} SRA experiment IDs, {len(projects)} BioProjects, "
          f"{len(ids)} DataCite DOI records, {len(titles)} exact titles and five source accessions.")


if __name__=="__main__":
    main()
