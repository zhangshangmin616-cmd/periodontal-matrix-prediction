"""Test the composition-only side contrast after equalizing day-7 UMI depth."""
import argparse
from pathlib import Path

import numpy as np

from benchmark import HERE, WORK, rows, write


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, default=WORK/'stage10/reference_counts.npz')
    parser.add_argument('--features', type=Path, default=WORK/'stage5/reference/features.tsv')
    parser.add_argument('--out', type=Path, default=HERE/'depth_sensitivity.tsv')
    a = parser.parse_args()
    features = rows(a.features)
    rng = np.random.default_rng(20261124)
    with np.load(a.counts, allow_pickle=False) as z:
        lookup = {g:i for i,g in enumerate(z['genes'].tolist())}
        ids = [lookup[r['gene']] for r in features]
        panel = np.array([r['ecm_component']=='other' for r in features])
        cp_profiles = z['CP'][ids] * (1e6 / z['CP_total'])
        predictions = []
        histograms = []
        for early, late in [('PP3D','PP7D'),('TP3D','TP7D')]:
            depth = z[early+'_total']
            predictions.append(np.log2(1 + cp_profiles @ (depth/depth.sum())))
            counts = z[late][ids].sum(axis=1)
            # Pool every other feature, preserving the true finite library size.
            histograms.append(np.r_[counts,z[late+'_total'].sum()-counts.sum()])
    p = predictions[0]-predictions[1]
    smaller = min(int(h.sum()) for h in histograms)
    results = []
    for depth in [smaller, smaller//2]:
        gains = {'background':[],'panel':[]}
        for _ in range(200):
            ys = []
            for h in histograms:
                draw = h if depth==h.sum() else rng.multivariate_hypergeometric(h,depth)
                ys.append(np.log2(1 + 1e6*draw[:-1]/depth))
            y = ys[0]-ys[1]
            gain = np.abs(y)-np.abs(y-p)
            gains['panel'].append(gain[panel].mean())
            gains['background'].append(gain.mean())
        for group, gs in gains.items():
            q = np.quantile(gs,[.025,.5,.975])
            results.append(dict(equal_umi_depth=depth, gene_set=group, draws=len(gs),
                q025_gain=q[0], median_gain=q[1], q975_gain=q[2],
                composition_beats_zero=int(np.sum(np.asarray(gs)>0)), seed=20261124))
    write(a.out, results)


if __name__=='__main__':
    main()
