"""Recover both Osaka arrays by matching scanner feature IDs to the platform."""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
EXT = WORK/'stage13/external'
sys.path.insert(0,str(EXT))
from extract import identifier_map
from extract_chip import soft

ACCESSIONS = ['GSM474118','GSM474198']


def main():
    raw=EXT/'raw'
    checks=json.loads((EXT/'chip_qc.json').read_text(encoding='utf8'))
    for name in ['GPL1291',*ACCESSIONS,*[s+'.txt.gz' for s in ACCESSIONS]]:
        assert hashlib.sha256((raw/name).read_bytes()).hexdigest()==checks[name+'_sha256']
    platform=soft('GPL1291','platform')
    ids=[r['ID'] for r in platform]
    idset=set(ids)
    assert len(ids)==len(idset)==30336
    with gzip.open(WORK/'stage5/paired/hgnc_ids.tsv.gz','rt',encoding='utf8') as f:
        entrez=identifier_map(csv.DictReader(f,delimiter='\t'),'entrez_id')
    scans=[]
    deposited=[]
    for gsm in ACCESSIONS:
        deposited.append({r['ID_REF']:float(r['VALUE']) for r in soft(gsm,'sample')})
        keyed={}
        with gzip.open(raw/(gsm+'.txt.gz'),'rt',encoding='utf8',errors='replace') as f:
            reader=csv.DictReader(f,delimiter='\t')
            reader.fieldnames=[k.strip() for k in reader.fieldnames]
            for r in reader:
                pid=r['ID']
                if pid not in idset:
                    continue
                assert pid not in keyed
                keyed[pid]=[float(r['Ch2 N Log Ratio']),float(r['Ch1 N (Median-B)']),float(r['Ch2 N (Median-B)'])]
        assert set(keyed)==idset
        scans.append(keyed)
    with gzip.open(HERE/'probes.tsv.gz','wt',encoding='utf8',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['probe','platform_symbol','symbol_by_entrez','logA','ch1A','ch2A','depositedA',
                    'logB','ch1B','ch2B','depositedB'])
        for r in platform:
            pid=r['ID']
            w.writerow([pid,r['Gene Symbol'] if 'Gene Symbol' in r else r['Symbol'],entrez.get(r['Entrez Gene ID'],''),
                        *scans[0][pid],deposited[0][pid],*scans[1][pid],deposited[1][pid]])


if __name__=='__main__':
    main()
