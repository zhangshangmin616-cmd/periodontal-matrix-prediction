"""Audit probe keys and retain measurements from the reconciled 48-hour array."""
import csv
import gzip
import hashlib
import json
from collections import defaultdict

import numpy as np

from extract import HERE, RAW, WORK, identifier_map
from score import save


def soft(s, kind):
    lines = (RAW/s).read_text(encoding="utf8").splitlines()
    i = lines.index("!"+kind+"_table_begin")
    j = lines.index("!"+kind+"_table_end")
    return list(csv.DictReader(lines[i+1:j], delimiter="\t"))


def extract():
    with gzip.open(WORK/"stage5/paired/hgnc_ids.tsv.gz","rt",encoding="utf8") as stream:
        mapping = identifier_map(csv.DictReader(stream,delimiter="\t"),"entrez_id")
    with (WORK/"stage5/reference/features.tsv").open(encoding="utf8") as stream:
        features = {r["gene"]: r["ecm_component"] for r in csv.DictReader(stream,delimiter="\t")}
    probes = {r["ID"]: mapping.get(r["Entrez Gene ID"]) for r in soft("GPL1291","platform")}
    audits, values = {}, []
    for gsm in ("GSM474118","GSM474198"):
        observations = {r["ID_REF"]: float(r["VALUE"]) for r in soft(gsm,"sample")}
        assert set(probes) == set(observations)
        with gzip.open(RAW/(gsm+".txt.gz"),"rt",encoding="utf8",errors="replace") as stream:
            reader = csv.DictReader(stream,delimiter="\t")
            reader.fieldnames = [f.strip() for f in reader.fieldnames]
            scan = {r["ID"]: (float(r["Ch2 N Log Ratio"]),float(r["Ch1 N (Median-B)"]),float(r["Ch2 N (Median-B)"]))
                    for r in reader if r["ID"] in probes}
        assert set(observations) == set(scan)
        exact = sum(v == scan[k][0] for k,v in observations.items())
        audits[gsm] = {"probes":len(scan),"processed_matches_attached_by_id":exact}
        if exact < len(scan):
            continue
        assert gsm == "GSM474198"
        record = (RAW/gsm).read_text(encoding="utf8")
        assert "#VALUE = normalized log2 ratio (Cy3/Cy5)" in record
        assert "!Sample_source_name_ch1 = 48-h stretched human periodontal ligament cell" in record
        assert "!Sample_label_ch1 = Cy3" in record
        for name, check_intensity in (("positive_signal", True),("all_published_values", False)):
            by_gene = defaultdict(list)
            for k, y in observations.items():
                gene = probes[k]
                if gene in features and (not check_intensity or min(scan[k][1:]) > 0):
                    by_gene[gene].append(y)
            for gene, vals in sorted(by_gene.items()):
                values.append(dict(rule=name,gene=gene,target=int(features[gene]=="other"),
                                   log2_change=np.mean(vals),probes=len(vals)))
    assert audits["GSM474118"]["processed_matches_attached_by_id"]<1000
    assert audits["GSM474198"]["processed_matches_attached_by_id"]==30336
    save("chip_effects.tsv",values)
    for name in ("GPL1291","GSM474118","GSM474198","GSM474118.txt.gz","GSM474198.txt.gz"):
        audits[name+"_sha256"] = hashlib.sha256((RAW/name).read_bytes()).hexdigest()
    (HERE/"chip_qc.json").write_text(json.dumps(audits,indent=2)+"\n",encoding="utf8")


if __name__ == "__main__":
    extract()
