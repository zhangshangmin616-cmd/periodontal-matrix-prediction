"""Independent Zhao patient-group forecasts and a separately invoked outcome evaluation."""
import argparse
import json
from datetime import datetime, timezone
import numpy as np
import pandas as pd
from forecast import CLASSES,ROOT,PAIR,OUT,Reference,filehash,put,score

ZHAO=ROOT/"zhao"
MODELS=["composition","state","zero","state_only"]

def layout(ref):
    runs=pd.read_csv(ZHAO/"joined_manifest.tsv",sep="\t").set_index("run")
    assert len(runs)==63 and runs.index.is_unique and runs["sample"].is_unique
    sizes=runs.groupby(["batch","group"]).size().to_dict()
    assert sizes=={(3,"control"):22,(3,"day28"):17,(4,"control"):5,(4,"day7"):19}
    ids=pd.read_csv(ZHAO/"feature_map.tsv",sep="\t",dtype={"entrez_id":str}).set_index("gene").loc[ref.genes]
    assert ids.entrez_id.is_unique and ids.in_all_runs.all()
    return runs,ids.entrez_id

def one_run(run,ids):
    df=pd.read_csv(ZHAO/"counts"/(run+".tsv.gz"),sep="\t",dtype={"Geneid":str}).set_index("Geneid")
    assert df.shape[1]==1 and df.index.is_unique
    v=df.iloc[:,0].astype(float)
    assert np.isfinite(v).all() and (v>=0).all() and v.sum()>0
    return v.loc[ids].to_numpy()*1e6/v.sum()

def aggregate(frame,cols,standardized):
    if standardized:
        assert set(frame.arch)=={"maxillary","mandibular"}
        return frame.groupby("arch")[cols].mean().mean()
    return frame[cols].mean()

def prepare():
    ref=Reference();runs,ids=layout(ref)
    control=runs[runs.group.eq("control")]
    count={r:one_run(r,ids) for r in control.index}
    forecast=[];w=[]
    for batch,group in control.groupby("batch"):
        for run,row in group.iterrows():
            x=count[run];peer=np.mean([count[s] for s in group.index if s!=run],axis=0)
            q,fit=ref.coefficients(x,peer)
            c,m=ref.predict(x,q,"primary")
            _,within=ref.predict(x,q,"state_only")
            forecast.append(dict(run=run,batch=batch,arch=row.arch,baseline=score(x,ref.idx),
                                 composition=c,state=m,zero=0.,state_only=within))
            w.extend(dict(run=run,cell_class=k,weight=v,anchor_rmse=fit) for k,v in zip(CLASSES,q))
    put(forecast,"external_control_forecasts.tsv");put(w,"external_weights.tsv")
    pred=[]
    for batch,frame in pd.DataFrame(forecast).groupby("batch"):
        for mode in ["batch","equal_jaw"]:
            shifts=aggregate(frame,MODELS,mode=="equal_jaw")
            pred.append(dict(batch=batch,time="day7" if batch==4 else "day28",analysis=mode,**shifts.to_dict()))
    put(pred,"external_predictions.tsv")
    locked=[ROOT/"external.py",ROOT/"external_plan.txt",ROOT/"external_secondary.txt",ROOT/"forecast.py",OUT/"lock.json",
        OUT/"external_control_forecasts.tsv",OUT/"external_weights.tsv",OUT/"external_predictions.tsv",ZHAO/"joined_manifest.tsv",
        ZHAO/"feature_map.tsv",ZHAO/"ncbi_run_manifest.tsv",*(ZHAO/"counts").glob("SRR*.tsv.gz")]
    (OUT/"external_lock.json").write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        inputs_before_scoring="Expression loaded for 27 controls only. Treated file bytes checked without parsing values.",
        sha256={p.relative_to(ROOT).as_posix():filehash(p) for p in locked}),indent=2))
    print(pd.DataFrame(pred).round(4).to_string(index=False))

def evaluate():
    from evaluate import verify
    verify()
    lock=json.loads((OUT/"external_lock.json").read_text())
    for f,sha in lock["sha256"].items():
        if filehash(ROOT/f)!=sha:raise ValueError(f"External lock mismatch for {f}")
    ref=Reference();runs,ids=layout(ref)
    for run in runs.index:
        runs.loc[run,"score"]=score(one_run(run,ids),ref.idx)
    runs.reset_index()[["run","group","batch","arch","sample","score"]].to_csv(OUT/"external_sample_scores.tsv",sep="\t",index=False)
    preds=pd.read_csv(OUT/"external_predictions.tsv",sep="\t")
    output=[]
    for row in preds.itertuples():
        ctrl=runs[(runs.batch==row.batch)&(runs.group=="control")]
        treated=runs[(runs.batch==row.batch)&(runs.group==row.time)]
        std=row.analysis=="equal_jaw"
        base=aggregate(ctrl,["score"],std)["score"];post=aggregate(treated,["score"],std)["score"]
        for key in MODELS:
            predicted=getattr(row,key)
            output.append(dict(time=row.time,batch=row.batch,analysis=row.analysis,model=key,
                n_control=len(ctrl),n_treated=len(treated),control_score=base,loaded_score=post,predicted=predicted,
                observed=post-base,abs_error=abs(predicted-(post-base))))
    put(output,"external_errors.tsv")
    controls=runs[runs.group=="control"]
    unmatched=[]
    for time in ["day7","day28"]:
        a=runs.loc[runs.group==time,"score"].mean();b=controls.score.mean()
        unmatched.append(dict(time=time,n_control=len(controls),loaded_mean=a,pooled_control_mean=b,
                              observed_unmatched=a-b))
    put(unmatched,"external_unmatched_controls.tsv")
    print(pd.DataFrame(output).query("analysis=='batch'")[["time","model","observed","predicted","abs_error"]].round(4).to_string(index=False))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("step",choices=["freeze","evaluate"])
    if p.parse_args().step=="freeze":prepare()
    else:evaluate()
