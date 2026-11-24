"""Build mapped culture count inputs from the eight checked SRA count files."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

import numpy as np


HERE=Path(__file__).resolve().parent
WORK=HERE.parents[1]


def main():
    with (WORK/'stage5/zhao/feature_map.tsv').open(encoding='utf8') as f:
        ids={r['gene']:r['entrez_id'] for r in csv.DictReader(f,delimiter='\t')}
    with np.load(WORK/'stage13/external/paired_counts.npz') as hk:
        names=hk['genes']
        target=hk['target']
    with (HERE/'count_manifest.tsv').open(encoding='utf8') as f:
        check={row['run']:row['md5'] for row in csv.DictReader(f,delimiter='\t')}
    info=json.loads((HERE/'ena_runs.json').read_text(encoding='utf8'))
    ids_order=[ids[name] for name in names]
    assert len(set(ids_order))==len(ids_order)
    profiles,totals,arms,runs=[],[],[],[]
    for row in sorted(info,key=lambda r:r['run_accession']):
        run=row['run_accession']
        with gzip.open(HERE/'counts'/(run+'.tsv.gz'),'rb') as f:
            raw=f.read()
        assert hashlib.md5(raw).hexdigest()==check[run]
        it=iter(raw.decode().splitlines())
        assert next(it).startswith('Geneid\t')
        entire=dict((a,float(b)) for a,b in (line.split('\t') for line in it))
        assert set(ids_order)<=set(entire)
        counts=np.array([entire[q] for q in ids_order])
        profiles.append(counts)
        totals.append(sum(entire.values()))
        freq=re.search(r'under (0(?:\.[15])?) Hz:',row['experiment_title']).group(1)
        arms.append(freq)
        runs.append(run)
    assert sorted(arms)==['0']*3+['0.1']*2+['0.5']*3
    np.savez_compressed(HERE/'xian_counts.npz',genes=names,target=target,counts=np.stack(profiles,axis=1),
                         total=np.array(totals),run=np.array(runs),arm=np.array(arms))


if __name__=='__main__':
    main()
