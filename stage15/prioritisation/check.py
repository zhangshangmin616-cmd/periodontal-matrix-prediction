"""Independently enumerate assignments and verify the 101-gene retrieval scores."""
import math
from shortlists import HERE, clinical, cultures, rows, top_weights
import numpy as np


def brute_weights(values, k):
    absolute = abs(values)
    result = []
    for a in absolute:
        larger = sum(1 for b in absolute if b > a)
        equals = sum(1 for b in absolute if b == a)
        result.append(min(1.,max(0.,(k-larger)/equals)))
    return result


def main():
    lookup = {(r["destination"],r["unit"],r["forecast"],int(r["k"])):r for r in rows(HERE/"scores.tsv")}
    verified = 0
    for dest, unit, model, genes, base, measured, forecast in [*clinical(),*cultures()]:
        cutpoints = np.percentile(sorted(base), [20,40,60,80])
        strata = [sum(v > c for c in cutpoints) for v in base]
        for k in (5,10,20):
            orig = lookup[dest,unit,model,k]
            p, y = brute_weights(forecast,k), brute_weights(measured,k)
            assert math.isclose(sum(p),k) and math.isclose(sum(y),k)
            count = sum(a*b for a,b in zip(p,y))
            signed = sum(a*b for i,(a,b) in enumerate(zip(p,y)) if measured[i]*forecast[i]>0)
            expected = expected_signed = 0.
            for i,wi in enumerate(p):
                members = [j for j in range(len(y)) if strata[i]==strata[j]]
                for j in members:
                    h = wi*y[j]/len(members)
                    expected += h
                    expected_signed += h * (forecast[i]*measured[j]>0)
            for key,num in zip(("overlap","signed_overlap","expected","expected_signed"),(count,signed,expected,expected_signed)):
                assert math.isclose(num,float(orig[key]),abs_tol=1.e-10), (unit, model, key)
            verified += 1
    assert verified == len(lookup)
    np.testing.assert_allclose(top_weights(np.array([2.,1.,-1.,1.,0.]),2), [1.,1/3,1/3,1/3,0.])
    print(f"Verified {verified} retrieval comparisons using all within-stratum pairings.")


if __name__=="__main__":
    main()
