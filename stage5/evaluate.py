"""Open exposed tissue counts only after checking the frozen forecast."""
import json
import numpy as np
import pandas as pd
from forecast import OUT, ROOT, Reference, filehash, load_bulk, meta, put, score

def verify():
    lock=json.loads((OUT/"lock.json").read_text())
    assert {"results/predictions.tsv","results/inner_predictions.tsv","design.txt"}<=set(lock["sha256"])
    for f,sha in lock["sha256"].items():
        if filehash(ROOT/f)!=sha: raise ValueError(f"Frozen file changed: {f}")

def main():
    verify()
    ref=Reference();info=meta()
    matrix=load_bulk(info.index.tolist(),ref.genes)
    info["score"]=[score(matrix[s].to_numpy(),ref.idx) for s in info.index]
    arms=info.groupby(["patient","arm"]).score.mean().unstack()
    observed=arms.force-arms.control
    single=info.groupby(["patient","arm"]).size().unstack().eq(1).all(axis=1)
    allpred=pd.read_csv(OUT/"predictions.tsv",sep="\t")
    primary=allpred.query("setting=='primary'").set_index("patient")
    inner=pd.read_csv(OUT/"inner_predictions.tsv",sep="\t")
    summary=[];errors=[];alpha=[]
    for patient,y in observed.items():
        train=inner[inner.held_out.eq(patient)].set_index("patient")
        assert patient not in train.index and set(train.index)==set(observed.index)-{patient}
        forecasts={"composition":primary.at[patient,"composition"],"state":primary.at[patient,"state"],
                   "zero":0.,"median_train":float(observed.drop(patient).median())}
        for model in ["composition","state"]:
            levels=np.arange(11)/10
            losses=np.abs(train[model].to_numpy()[:,None]*levels-observed.loc[train.index].to_numpy()[:,None]).mean(axis=0)
            factor=levels[np.nonzero(losses<=losses.min()+1e-12)[0][0]]
            forecasts["scaled_"+model]=factor*forecasts[model]
            alpha.append(dict(patient=patient,model=model,multiplier=factor,training_mae=losses.min()))
        for model,p in forecasts.items():
            errors.append(dict(patient=patient,model=model,predicted=p,observed=y,abs_error=abs(p-y),
                single_tooth_each=bool(single.loc[patient]),arch_tooth_matched=patient in [3,7,8,9,10]))
    e=pd.DataFrame(errors)
    for subset,df in [("all",e),("single_tooth_each",e[e.single_tooth_each]),("arch_tooth_matched",e[e.arch_tooth_matched])]:
        for name,g in df.groupby("model",sort=False):
            summary.append(dict(subset=subset,model=name,n=len(g),mae=g.abs_error.mean(),
                median_absolute_error=g.abs_error.median(),median_prediction=g.predicted.median(),mean_prediction=g.predicted.mean(),
                median_observed=g.observed.median(),mean_observed=g.observed.mean(),mean_bias=(g.predicted-g.observed).mean(),
                direction_agrees=int((np.sign(g.predicted)==np.sign(g.observed)).sum())))
    sens=[];pairs=[]
    for setting,rows in allpred.groupby("setting",sort=False):
        p=rows.set_index("patient").loc[observed.index]
        c=(p.composition-observed).abs();m=(p.state-observed).abs()
        for name in ["composition","state"]:
            dif=p[name]-observed
            sens.append(dict(setting=setting,model=name,n=len(p),mae=dif.abs().mean(),mean_bias=dif.mean(),
                             zero_mae=observed.abs().mean()))
        for i in p.index:
            pairs.append(dict(setting=setting,patient=i,composition_error=c[i],state_error=m[i],
                              gain_over_composition=c[i]-m[i],state_gain_over_zero=abs(observed[i])-m[i]))
    put(errors,"patient_errors.tsv");put(summary,"paired_metrics.tsv");put(alpha,"calibration.tsv")
    put(sens,"settings_metrics.tsv");put(pairs,"patient_gains.tsv")
    info.reset_index().to_csv(OUT/"observed_teeth.tsv",sep="\t",index=False)
    print(pd.DataFrame(summary).query("subset=='all'")[["model","n","mae","mean_prediction","mean_observed"]].round(4).to_string(index=False))

if __name__=="__main__":main()
