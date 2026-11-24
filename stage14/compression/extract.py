"""Rebuild compact static-compression count matrices from public count files."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import tarfile
from collections import Counter
from io import TextIOWrapper
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
URLS = {
    'GSE112122_RAW.tar': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112122/suppl/GSE112122_RAW.tar',
    'GSE109167_raw_counts_GRCh38.p13_NCBI.tsv.gz': 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109167&format=file&type=rnaseq_counts&file=GSE109167_raw_counts_GRCh38.p13_NCBI.tsv.gz',
    'GSE347200_gene_count_matrix.txt.gz': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE347nnn/GSE347200/suppl/GSE347200_gene_count_matrix.txt.gz',
}


def registry_map(path, field):
    with gzip.open(path,'rt',encoding='utf8') as handle:
        pairs = {(r['symbol'],r[field]) for r in csv.DictReader(handle,delimiter='\t') if r[field]}
    symbols, ids = map(Counter,zip(*pairs))
    return {k:v for k,v in pairs if symbols[k] == ids[v] == 1}


def matrix_table(path, samples, keys, ensembl):
    with gzip.open(path,'rt') as handle:
        stream = csv.reader(handle,delimiter='\t')
        header = next(stream)
        columns = [header.index(s) for s in samples]
        index = {}
        for row in stream:
            gene = row[0].split('.')[0] if ensembl else row[0]
            if ensembl and not gene.startswith('ENSG'):
                continue
            assert gene not in index
            index[gene] = [float(row[c]) for c in columns]
    return np.asarray([index[k] for k in keys]).T, np.sum(list(index.values()),axis=0)


def bangkok(path, keys):
    profiles, totals = {}, {}
    with tarfile.open(path) as archive:
        for member in archive:
            m = re.fullmatch(r'GSM\d+_count([123])(St|Cs)\.txt\.gz',member.name)
            if m is None:
                continue
            with TextIOWrapper(gzip.GzipFile(fileobj=archive.extractfile(member))) as handle:
                entries = [(r[0].split('.')[0],int(r[1])) for r in csv.reader(handle,delimiter='\t') if r[0].startswith('ENSG')]
            index = dict(entries)
            assert len(index) == len(entries)
            label = m[2]+m[1]
            profiles[label] = [index[k] for k in keys]
            totals[label] = sum(index.values())
    names = ['Cs1','Cs2','Cs3','St1','St2','St3']
    return names, np.asarray([profiles[k] for k in names]), np.array([totals[k] for k in names])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir',type=Path,required=True)
    args = parser.parse_args()
    features = list(csv.DictReader((WORK/'stage5/reference/features.tsv').read_text().splitlines(),delimiter='\t'))
    genes = np.array([r['gene'] for r in features])
    targets = np.array([r['ecm_component']=='other' for r in features])
    assert len(genes)==10806 and targets.sum()==101
    registry = WORK/'stage5/paired/hgnc_ids.tsv.gz'
    ensg = registry_map(registry,'ensembl_gene_id')
    entrez = registry_map(registry,'entrez_id')
    eid = [ensg[g] for g in genes]
    nid = [entrez[g] for g in genes]
    filenames = list(URLS)
    bang_samples,bang_counts,bang_totals = bangkok(args.raw_dir/filenames[0],eid)
    beijing_samples = ['GSM2934477','GSM2934478','GSM2934479','GSM2934480']
    bj_counts,bj_totals = matrix_table(args.raw_dir/filenames[1],beijing_samples,nid,False)
    xian_samples = ['C_1','C_2','C_3','A_1','A_2','A_3']
    xa_counts,xa_totals = matrix_table(args.raw_dir/filenames[2],xian_samples,eid,True)
    counts = np.vstack([bang_counts,bj_counts,xa_counts])
    totals = np.concatenate([bang_totals,bj_totals,xa_totals])
    assert (counts>=0).all() and (counts.sum(axis=1)<=totals).all()
    np.savez_compressed(HERE/'counts.npz',gene=genes,target=targets,
                        accession=np.repeat(['GSE112122','GSE109167','GSE347200'],[6,4,6]),
                        arm=np.array(['control']*3+['force']*3+['control']*2+['force']*2+['control']*3+['force']*3),
                        sample=np.array(bang_samples+beijing_samples+xian_samples),
                        count=counts,library_total=totals)
    hashes = {n:{'url':URLS[n],'sha256':hashlib.sha256((args.raw_dir/n).read_bytes()).hexdigest()} for n in filenames}
    (HERE/'sources.json').write_text(json.dumps(hashes,indent=2)+'\n',encoding='utf8')


if __name__=='__main__':
    main()
