"""Score both RNA assays from Masson (2013) on the frozen day-seven forecasts."""
import csv
import gzip
import json
import math
from pathlib import Path
from statistics import fmean


HERE = Path(__file__).resolve().parent
SOURCE = 'https://repositorio.pgsscogna.com.br/bitstream/123456789/2736/1/282d74066b1278126458adace089b49c.pdf'
REGIONS = ['VC', 'VM', 'VA', 'LA', 'LM', 'LC']


def read_source_table(name):
    with (HERE / 'human' / name).open(encoding='utf8') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_table(name, rows):
    with (HERE / name).open('w', encoding='utf8', newline='') as stream:
        writer = csv.DictWriter(stream, rows[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    with gzip.open(HERE / 'state_multipliers.tsv.gz', 'rt') as stream:
        forecasts = {r['gene']: r for r in csv.DictReader(stream, delimiter='\t')
                     if r['gene'] in {'IFNG', 'IL10', 'ACTB'}}
    sites = []
    figures = []
    regions = {(r['gene'], r['anatomical_site'], int(r['day'])): float(r['median'])
               for r in read_source_table('masson_anatomical_medians.tsv')}
    groups = read_source_table('masson_zone_day7.tsv')
    assert {r['gene'] for r in groups} == {'IFNG', 'IL10'}
    for row in groups:
        gene = row['gene']
        pressure, tension = [float(forecasts[gene][side]) - float(forecasts['ACTB'][side])
                             for side in ['PP7D', 'TP7D']]
        contrast = pressure - tension
        c, t, middle = (float(row[x]) for x in ['pressure_mean', 'tension_mean', 'intermediate_mean'])
        for scale, y in [('linear_relative_expression', math.log2(c / t)),
                         ('already_log2_expression_sensitivity', c - t)]:
            figures.append(dict(gene=gene, scale=scale, pressure_mean=c, tension_mean=t,
                                intermediate_mean=middle, predicted_C_minus_T=contrast,
                                observed_C_minus_T=y, abs_model_error=abs(y-contrast),
                                abs_zero_error=abs(y), gain=abs(y) - abs(y-contrast)))
        for site in REGIONS:
            c0, treated = (regions[gene, site, d] for d in (0, 7))
            assert treated > c0
            sites.append(dict(gene=gene, region=site, untreated_median=c0, day7_median=treated,
                              log2_ratio_if_linear=math.log2(treated / c0),
                              forecast_if_pressure=pressure, forecast_if_tension=tension,
                              both_forecast_directions_opposite=pressure < 0 and tension < 0))
    write_table('masson_sides.tsv', figures)
    write_table('masson_sites.tsv', sites)
    scales = {}
    for scale in sorted({r['scale'] for r in figures}):
        subset = [r for r in figures if r['scale'] == scale]
        scales[scale] = {name: fmean(r[name] for r in subset)
                         for name in ['abs_model_error', 'abs_zero_error', 'gain']}
    data = dict(source=SOURCE, n_teeth_day7=14, n_teeth_control=20,
                n_patients_all_timepoints=33, n_teeth_all_timepoints=84,
                results_by_assumed_scale=scales)
    (HERE / 'masson_summary.json').write_text(json.dumps(data, indent=2) + '\n', encoding='utf8')
    print(scales)


if __name__ == '__main__':
    main()
