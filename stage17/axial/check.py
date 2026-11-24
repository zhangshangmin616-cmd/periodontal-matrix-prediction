"""Independently verify group-effect and gene-error summaries of the Chengdu contrasts."""
import csv
import gzip
import json
from math import fsum, isclose
from pathlib import Path

HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def main():
    with gzip.open(WORK/'stage8/axial_log2_multipliers.tsv.gz','rt') as stream:
        preds={r['gene']:r for r in csv.DictReader(stream,delimiter='\t') if r['in_panel']=='True'}
    with (HERE/'scores.tsv').open(encoding='utf8') as f:
        output={r['condition']:r for r in csv.DictReader(f,delimiter='\t')}
    assert len(preds)==101 and set(output)=={'I3','I10','E3','E10'}
    for name,result in output.items():
        path=WORK/f'stage8/hotm/deg_C0_{name}.json'
        doc=json.loads(path.read_bytes())
        metadata=json.loads(path.with_name(path.stem+'_meta.json').read_bytes())
        assert f'Select Group 1:\nC0\nSelect Group 2:\n{name}\n' in metadata['label']
        effects={r[1]:r[3] for r in doc['rows'] if r[1] in preds and r[3] is not None}
        n=len(effects)
        residuals=[float(preds[g][name])-v for g,v in effects.items()]
        obs_sum=fsum(effects.values())
        panel_err=abs(fsum(residuals))/n
        zero_panel_err=abs(obs_sum)/n
        gene_err=fsum(map(abs,residuals))/n
        zero_gene_err=fsum(map(abs,effects.values()))/n
        for field,reference in [('measured_zero_panel_error',zero_panel_err),('measured_model_panel_error',panel_err),
                                ('measured_zero_gene_mae',zero_gene_err),('measured_model_gene_mae',gene_err)]:
            assert isclose(float(result[field]),reference,abs_tol=1e-12),name
        lo,hi=float(result['full_panel_gain_lo']),float(result['full_panel_gain_hi'])
        if n!=101:
            assert set(preds)-set(effects)=={'SGCZ'}
            p=fsum(float(r[name]) for r in preds.values())
            extrema=[(abs(obs_sum+y)-abs(p-obs_sum-y))/101 for y in (1e3,-1e3)]
            assert isclose(lo,min(extrema),abs_tol=1e-12) and isclose(hi,max(extrema),abs_tol=1e-12)
        else:
            assert isclose(hi,zero_panel_err-panel_err,abs_tol=1e-12) and lo==hi
    print("Verified four Chengdu contrasts against original metadata and unthresholded coefficients.")


if __name__=='__main__':
    main()
