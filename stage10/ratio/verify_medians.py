"""Check two majority order-statistic bounds and the locked model error difference."""
import argparse
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen


HERE=Path(__file__).resolve().parent
FIGURE='https://cdn.ncbi.nlm.nih.gov/pmc/blobs/7c18/2719712/c6654d00da6e/nihms125069f3.jpg'
DIGEST='120985cab274a28ace6f36ea8fde032340071a3e4facb49c62968c7472d0398c'


def read(path):
    with path.open(encoding='utf8') as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--predictions',type=Path,default=HERE.parents[1]/'stage6/target_predictions.tsv')
    parser.add_argument('--figure',type=Path,help='Optional saved Fig. 3 source file')
    parser.add_argument('--fetch-figure',action='store_true',help='Check the transcribed source image by SHA-256')
    parser.add_argument('--out',type=Path,default=HERE/'median_check.json')
    a=parser.parse_args()
    if a.figure or a.fetch_figure:
        if a.figure:
            img=a.figure.read_bytes()
        else:
            with urlopen(Request(FIGURE,headers={'User-Agent':'Mozilla/5.0'}),timeout=30) as r:
                img=r.read()
        assert hashlib.sha256(img).hexdigest()==DIGEST
    counts=read(HERE/'figure3_groups.tsv')
    for side,key in [('pressure','[1,+inf)'),('tension','(-inf,-1]')]:
        group=[r for r in counts if r['side']==side]
        assert {int(r['denominator']) for r in group}=={26}
        assert sum(int(r['cases']) for r in group)==26
        assert sum(int(r['cases']) for r in group if r['class_interval_log2_rankl_over_opg']==key)>=14
    pred={(r['side'],r['gene']):Decimal(r['predicted_log2_ratio']) for r in read(a.predictions)}
    pressure=pred['pressure','TNFSF11']-pred['pressure','TNFRSF11B']
    tension=pred['tension','TNFSF11']-pred['tension','TNFRSF11B']
    delta=pressure-tension
    bound=Decimal('2')
    assert 0<delta<bound
    for y in (bound,Decimal('3'),Decimal('100')):
        assert abs(y)-abs(y-delta)==delta
    result=dict(figure_url=FIGURE,figure_sha256=DIGEST,freeze_predictions_sha256=hashlib.sha256(a.predictions.read_bytes()).hexdigest(),
                observed_lower_bound=str(bound),pressure_prediction=str(pressure),tension_prediction=str(tension),
                predicted_median_contrast=str(delta),absolute_error_gain_vs_no_side_difference=str(delta),
                minimum_model_abs_error=str(bound-delta))
    a.out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':
    main()
