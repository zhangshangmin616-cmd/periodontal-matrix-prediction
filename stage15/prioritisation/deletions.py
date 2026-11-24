"""Score the same clinical shortlists after deleting one contributing observation."""
import sys
import importlib.util
import numpy as np
from shortlists import HERE, WORK, clinical, evaluate, write


def main():
    prior = WORK/"stage13/transfer/deletions.py"
    spec = importlib.util.spec_from_file_location("old_deletions", prior)
    old = importlib.util.module_from_spec(spec)
    sys.path.insert(0,str(prior.parent))
    spec.loader.exec_module(old)
    with np.load(WORK/"stage11/specificity/rank_changes.npz") as r:
        g, mask = r["gene"], r["component"]=="other"
        ranks = old.clinical_ranks(g,mask)
        full, drops = old.observed_leuven(ranks, "day7")
        order = np.argsort(g[mask])
        source = r["observed"][np.char.startswith(r["units"],"Kim_")][:,mask][:,order]
        units = r["units"][np.char.startswith(r["units"],"Kim_")]
    samples = []
    for dst, unit, model, genes, b, y, p in clinical():
        if dst != "Leuven_day7":
            continue
        assert np.array_equal(genes, g[mask][order])
        np.testing.assert_allclose(y,full[order],rtol=0,atol=1e-10)
        for role, run, observed in drops:
            if role=="Leuven_exposed_specimen":
                samples.append((dst,run,model,genes,b,observed[order],p))
        if model=="other_clinical":
            np.testing.assert_allclose(source.mean(axis=0),p,rtol=0,atol=1e-10)
            for i, donor in enumerate(units):
                samples.append((dst,donor,model,genes,b,y,np.delete(source,i,axis=0).mean(axis=0)))
    results,_ = evaluate(samples)
    assert len(results)==(2*19+11)*3
    write(HERE/"deletions.tsv",results)


if __name__=="__main__":
    main()
