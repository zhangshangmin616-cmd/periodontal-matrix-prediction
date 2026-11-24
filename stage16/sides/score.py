"""Check observational equivalence of three side assignments in the original clinical forecast."""
import csv
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def ranks(a):
    output=[]
    for vector in a:
        x=np.sort(vector)
        output.append(50*(x.searchsorted(vector)+x.searchsorted(vector,side="right")-1)/(len(x)-1))
    return np.array(output)


def write(path,records):
    with path.open("w",encoding="utf8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t",lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def main():
    with np.load(WORK/"stage8/mixture_inputs.npz",allow_pickle=False) as v:
        b,p,t,genes,ix=[v[k] for k in ("control","pressure","tension","symbols","idx")]
    assert p.shape==t.shape==b.shape==(5,10806) and len(ix)==101
    with np.load(WORK/"stage11/specificity/rank_changes.npz",allow_pickle=False) as v:
        unit=np.flatnonzero(v["units"]=="Zhao_day7").item()
        gene_index={g:i for i,g in enumerate(v["gene"])}
        sel=[gene_index[g] for g in genes[ix]]
        assert (v["component"][sel]=="other").all()
        measured=v["observed"][unit,sel]
        previous=v["fitted"][unit,sel]
    m=(p+t)*.5
    ref=ranks(b)[:,ix]
    metrics,by_gene=[],[]
    for model,(left,right) in {"original":(p,t),"side_labels_exchanged":(t,p),"same_on_both_sides":(m,m)}.items():
        pooled=(left+right)*.5
        assert np.array_equal(pooled,m)
        estimate=(ranks(pooled)[:,ix]-ref).mean(axis=0)
        np.testing.assert_allclose(estimate,previous,atol=2e-12,rtol=0)
        pressure=(ranks(left)[:,ix]-ref).mean(axis=0)
        tension=(ranks(right)[:,ix]-ref).mean(axis=0)
        contrast=pressure-tension
        delta=estimate-measured
        metrics.append(dict(assignment=model,targets=len(ix),unit="expression percentile points",
                            observed_panel_change=measured.mean(),predicted_panel_change=estimate.mean(),
                            panel_error=abs(delta.mean()),gene_mae=np.abs(delta).mean(),
                            zero_panel_error=abs(measured.mean()),zero_gene_mae=np.abs(measured).mean(),
                            predicted_side_mean_absolute_contrast=np.abs(contrast).mean(),
                            targets_with_side_contrast=int(np.sum(np.abs(contrast)>1e-10))))
        for i,g in enumerate(genes[ix]):
            by_gene.append(dict(assignment=model,gene=g,
                                predicted_pressure=pressure[i],predicted_tension=tension[i],
                                predicted_pressure_minus_tension=contrast[i],
                                predicted_pooled=estimate[i],observed_unsplit_change=measured[i]))
    write(HERE/"metrics.tsv",metrics)
    write(HERE/"forecasts.tsv",by_gene)
    print("Checked three assignments on",len(ix),"targets across",len(b),"clinical baselines.")


if __name__=="__main__":
    main()
