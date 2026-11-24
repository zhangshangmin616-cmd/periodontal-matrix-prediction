"""Learn a single temporal transfer coefficient in the source atlas, then carry an
external cohort's day-3 force contrast to day 10 without external calibration."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import numpy as np
from score_direct import metrics, write


HERE = Path(__file__).resolve().parent


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--atlas', type=Path, default=HERE.parents[1]/'stage8/axial_log2_multipliers.tsv.gz')
    parser.add_argument('--fit-only', action='store_true')
    args = parser.parse_args()
    with gzip.open(args.atlas, 'rt', encoding='utf8') as stream:
        atlas = list(csv.DictReader(stream, delimiter='\t'))
    training = [r for r in atlas if r['in_panel']=='False']
    train_x = np.array([float(r['I3'])-float(r['E3']) for r in training])
    train_y = np.array([float(r['I10'])-float(r['E10']) for r in training])
    # I10/E10 in the historical external-transfer file denote the day-7 atlas state.
    # Coefficient is unconstrained with no intercept. No external study data needed.
    slope = float(train_x@train_y / (train_x@train_x))
    fit = dict(training_genes=len(training), coefficient=slope, model='ordinary least squares, no intercept',
               atlas_days=[3,7], external_days=[3,10], external_outcomes_used_for_fit=False,
               atlas_sha256=hashlib.sha256(args.atlas.read_bytes()).hexdigest())
    (HERE / 'atlas_persistence_fit.json').write_text(json.dumps(fit, indent=2)+'\n', encoding='utf8')
    print(fit)
    if args.fit_only:
        return
    with gzip.open(HERE/'direct_log2fc.tsv.gz', 'rt', encoding='utf8') as stream:
        observations = [r for r in csv.DictReader(stream, delimiter='\t') if r['day3'] and r['day10']]
    observed = {r['gene']:float(r['day10']) for r in observations}
    day3 = {r['gene']:float(r['day3']) for r in observations}
    panel = [r['gene'] for r in atlas if r['in_panel']=='True' and r['gene'] in observed]
    rows = []
    scales = [('zero',0.), ('persistence',1.), ('atlas_scaled_persistence',slope), ('quarter_persistence',.25), ('half_persistence',.5)]
    for group, genes in [('fixed_background', list(observed)), ('fixed_panel',panel)]:
        for method, scale in scales:
            forecast = {g:scale*day3[g] for g in genes}
            rows.append(metrics(genes, observed, forecast, method, group,10))
    write(HERE/'atlas_persistence_metrics.tsv', rows)
    print([(r['group'],r['method'],r['mae'],r['gain_vs_zero']) for r in rows])
    # Printed means from the separately recruited Brazilian study, Table 3 on PDF p32.
    # Do not add the thesis' two measures to the transcriptome-panel accuracy.
    clinical = [('IFNG',5.50,3.26,7.24,4.08), ('IL10',4.53,3.08,6.37,3.92)]
    per_gene = []
    for label, transform in [('linear_relative', lambda p,t: np.log2(p/t)), ('already_log2', lambda p,t: p-t)]:
        for gene,p3,t3,p7,t7 in clinical:
            x,y = transform(p3,t3), transform(p7,t7)
            for method, scale in scales:
                forecast = scale*x
                per_gene.append(dict(scale=label,gene=gene,day3_pressure=p3,day3_tension=t3,
                                     day7_pressure=p7,day7_tension=t7,method=method,day3_observed=x,
                                     day7_observed=y,forecast=forecast,absolute_error=abs(y-forecast)))
    write(HERE/'masson_time_transfer.tsv', per_gene)
    for label in ['linear_relative','already_log2']:
        for method,_ in scales:
            errors = [r['absolute_error'] for r in per_gene if r['scale']==label and r['method']==method]
            print('Masson', label,method,np.mean(errors))


if __name__=='__main__':
    main()
