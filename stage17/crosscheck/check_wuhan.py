"""Recalculate all Wuhan endpoints from independent sample mapping and integer atlas counts."""
import csv
import gzip
import hashlib
import io
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
STAGE = HERE.parent
WORK = STAGE.parent
ORIG = STAGE / "wuhan_counts"


def read(path):
    with path.open(encoding="utf8",newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def midranks(profile):
    order = np.argsort(profile,kind="stable")
    x = profile[order]
    boundaries = np.r_[0, 1+np.flatnonzero(x[1:]!=x[:-1]), len(x)]
    mid = (boundaries[:-1]+boundaries[1:]-1)*(50.0/(len(x)-1))
    ranks = np.empty(len(x))
    ranks[order] = np.repeat(mid,np.diff(boundaries))
    return ranks


def close(a,b):
    np.testing.assert_allclose(a,b,atol=2e-10,rtol=0)


def predictors(symbols):
    with np.load(WORK/"stage10/reference_counts.npz",allow_pickle=False) as r:
        row = {s:i for i,s in enumerate(r["genes"])}
        sel = [row[s] for s in symbols]
        control = 1e6*r["CP"][sel]/r["CP_total"]
        loaded = 1e6*r["TP3D"][sel]/r["TP3D_total"]
        fibro = list(r["classes"]).index("fibroblast_like")
        residual = list(r["classes"]).index("residual")
        loaded[:,residual] = control[:,residual]
        w = r["CP_total"]/r["CP_total"].sum()
        f = np.log2((loaded[:,fibro]+1)/(control[:,fibro]+1))
        t = np.log2((1+(loaded*w).sum(axis=1))/(1+(control*w).sum(axis=1)))
    with gzip.open(WORK/"stage8/axial_log2_multipliers.tsv.gz","rt",encoding="utf8") as fp:
        old = {r["gene"]:float(r["E3"]) for r in csv.DictReader(fp,delimiter="\t")}
    close(t,[old[g] for g in symbols])
    return {"tissue":t,"fibroblast":f}


def samples(entrez):
    ena = json.loads((STAGE/"culture/ena_runs.json").read_bytes())
    mapping = defaultdict(list)
    hours = {}
    for run in ena:
        label = run["experiment_title"]
        m = re.search(r"(GSM[0-9]+): sh-Control, tension (0|24|48) h, rep [1-3]",label)
        if m:
            mapping[m[1]].append(run["run_accession"])
            hours[m[1]] = int(m[2])
    assert len(mapping)==9 and sum(map(len,mapping.values()))==12
    hashes = {r["run"]:r["md5"] for r in read(ORIG/"count_manifest.tsv")}
    decimal_counts, float_counts, totals = {}, {}, {}
    for sample,runs in mapping.items():
        cents = np.zeros(len(entrez),dtype=np.int64)
        counts = np.zeros(len(entrez))
        tot = 0
        for accession in sorted(runs):
            with gzip.open(ORIG/"counts"/(accession+".tsv.gz"),"rb") as archive:
                data = archive.read()
            assert hashlib.md5(data).hexdigest()==hashes[accession]
            table = csv.reader(io.StringIO(data.decode("utf8")),delimiter="\t")
            assert next(table)==["Geneid",accession+"_count"]
            raw = dict(table)
            def hundredths(value):
                whole, mark, frac = value.partition(".")
                assert not mark or len(frac)==2
                return 100*int(whole)+(int(frac) if mark else 0)
            exact = {g:hundredths(s) for g,s in raw.items()}
            cents += [exact[i] for i in entrez]
            counts += [float(raw[i]) for i in entrez]
            tot += sum(exact.values())
        totals[sample] = tot/100.
        decimal_counts[sample] = cents/100.
        float_counts[sample] = counts
    return hours,decimal_counts,float_counts,totals


def score(data,totals,times,coefs,target):
    cpm = {s:v*1e6/totals[s] for s,v in data.items()}
    ranks = {s:midranks(v)[target] for s,v in cpm.items()}
    logs = {s:np.log2(1+v[target]) for s,v in cpm.items()}
    rows,details = {},{}
    def average(arr,ids):
        return np.mean([arr[s] for s in ids],axis=0)
    for h in (24,48):
        eligible = sorted(s for s in cpm if times[s] in (0,h))
        for model,p in coefs.items():
            rank_predictions = {s:midranks(v*np.exp2(p))[target]-ranks[s] for s,v in cpm.items() if times[s]==0}
            for omit in ["",*eligible]:
                base = [s for s in eligible if times[s]==0 and s!=omit]
                loaded = [s for s in eligible if times[s]==h and s!=omit]
                for scale,values in [("rank",ranks),("log2",logs)]:
                    obs = average(values,loaded)-average(values,base)
                    pred = average(rank_predictions,base) if scale=="rank" else p[target]
                    modes = [(model,pred)]
                    if model=="tissue":
                        modes.append(("tissue_common_shift",np.full(len(obs),pred.mean())))
                    for mode,pp in modes:
                        key = (h,mode,scale,omit)
                        pe,ze = abs(pp.mean()-obs.mean()),abs(obs.mean())
                        ge,zg = np.mean(abs(pp-obs)),np.mean(abs(obs))
                        rows[key] = dict(observed_mean=obs.mean(),predicted_mean=pp.mean(),panel_abs_error=pe,
                                         zero_panel_abs_error=ze,panel_gain=ze-pe,gene_mae=ge,zero_gene_mae=zg,gene_gain=zg-ge)
                        if not omit:
                            details[key[:3]]=(obs,pp)
    return rows,details


def main():
    with np.load(WORK/"stage11/specificity/rank_changes.npz",allow_pickle=False) as z:
        genes, target = z["gene"],z["component"]=="other"
    ids = {r["gene"]:r["entrez_id"] for r in read(WORK/"stage5/zhao/feature_map.tsv")}
    times,exact,floats,totals = samples([ids[g] for g in genes])
    coefs = predictors(genes)
    rows,detail = score(floats,totals,times,coefs,target)
    compare,other = score(exact,totals,times,coefs,target)
    expected = read(ORIG/"scores.tsv")+read(ORIG/"deletions.tsv")
    assert len(expected)==len(rows)==84
    for r in expected:
        k=(int(r["hours"]),r["model"],r["scale"],r.get("omitted",""))
        for name,value in rows[k].items():
            close(value,float(r[name]))
    idx = {s:i for i,s in enumerate(genes[target])}
    detailed = read(ORIG/"targets.tsv")
    assert len(detailed)==12*101 and len(detail)==12
    for r in detailed:
        obs,p=detail[int(r["hours"]),r["model"],r["scale"]]
        j=idx[r["gene"]]
        close([obs[j],p[j]],[float(r["observed"]),float(r["predicted"])])
    summary = []
    for h in (24,48):
        k = (h,"tissue","rank","")
        d = [v for kk,v in rows.items() if kk[:3]==k[:3] and kk[3]]
        summary.append(dict(hours=h,zero_panel_error=rows[k]["zero_panel_abs_error"],model_panel_error=rows[k]["panel_abs_error"],
                            zero_gene_error=rows[k]["zero_gene_mae"],model_gene_error=rows[k]["gene_mae"],
                            deletions_with_panel_gain=sum(v["panel_gain"]>0 for v in d),
                            deletions_with_gene_loss=sum(v["gene_gain"]<0 for v in d),n_deletions=len(d),
                            exact_decimal_panel_gain=compare[k]["panel_gain"],exact_decimal_gene_gain=compare[k]["gene_gain"]))
    with (HERE/"wuhan_verified.tsv").open("w",encoding="utf8",newline="") as handle:
        w=csv.DictWriter(handle,fieldnames=list(summary[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(summary)
    print("Reproduced 84 endpoint/deletion rows and 1,212 gene-level records.")


if __name__=="__main__":
    main()
