"""Extract only the fixed panel's published, selection-conditioned culture readouts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET


HERE=Path(__file__).resolve().parent
W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
BASE='https://ars.els-cdn.com/content/image/1-s2.0-S0171933523000341-'


def data(path):
    with ZipFile(path) as z:
        document=ET.fromstring(z.read('word/document.xml'))
    def text(el):
        return ''.join(v.text or '' for v in el.iter(W+'t')).strip()
    table=document.find('.//'+W+'tbl')
    return [[text(c) for c in row.findall(W+'tc')] for row in table.findall(W+'tr')]


def main():
    cli=argparse.ArgumentParser()
    cli.add_argument('--source',type=Path,default=HERE)
    cli.add_argument('--features',type=Path,default=HERE.parents[1]/'stage5/reference/features.tsv')
    args=cli.parse_args()
    with args.features.open(encoding='utf8') as f:
        genes=[r['gene'] for r in csv.DictReader(f,delimiter='\t') if r['ecm_component']=='other']
    assert len(set(genes))==len(genes)==101
    a,b=[data(args.source/f'mmc{k}.docx') for k in [4,5]]
    assert a[0]==['Gene name','Gene biotype','log2-fold change','adjustedp-value']
    assert b[0][-6:]==[f'{g} {i}' for g in ['Control','Tension'] for i in range(1,4)]
    coef={r[0]:r[2] for r in a[1:]}
    counts={r[0]:r[-6:] for r in b[1:]}
    assert len(coef)==len(a)-1==1336 and len(counts)==len(b)-1==124
    assert counts.keys()<=coef.keys()
    with (HERE/'culture_readouts.tsv').open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['gene','reported_log2','control1','control2','control3','tension1','tension2','tension3'])
        for gene in genes:
            w.writerow([gene,coef.get(gene,''),*counts.get(gene,['']*6)])
    manifest={}
    for i in [4,5]:
        name=f'mmc{i}.docx'
        manifest[name]={'url':BASE+name,'sha256':hashlib.sha256((args.source/name).read_bytes()).hexdigest()}
    (HERE/'source.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    main()
