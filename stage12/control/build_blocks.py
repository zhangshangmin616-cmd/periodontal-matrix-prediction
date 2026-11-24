"""Rebuild disjoint cell-block counts from public GSE327337 matrices and saved gates."""
import argparse
import csv
import gzip
import hashlib
from pathlib import Path

import numpy as np
from scipy import io, sparse


HERE = Path(__file__).resolve().parent
WORK = HERE.parents[1]
LIBS = ['CP', 'PP3D', 'TP3D', 'PP7D', 'TP7D']
CLASSES = ['fibroblast_like', 'mural', 'immune', 'endothelial', 'glial', 'residual']


def read_rows(path, header=True):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t') if header else csv.reader(stream, delimiter='\t')
        return list(reader)


def write_table(path, records):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', encoding='utf8', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=records[0].keys(), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)


def materialize_membership(raw, flags):
    membership = []
    for lib in LIBS:
        original = read_rows(next(raw.glob(f'*_{lib}_barcodes.tsv.gz')), header=False)
        audit = read_rows(flags / f'{lib}_barcodes.tsv.gz')
        assert [r[0] for r in original] == [r['barcode'] for r in audit]
        accepted = [r for r in audit if r['passes_analysis'].lower() == 'true']
        assert len(set(r['barcode'] for r in accepted)) == len(accepted)
        for group in CLASSES:
            barcodes = [r['barcode'] for r in accepted
                        if (r['exclusive'] if r['exclusive'] in CLASSES[:-1] else 'residual') == group]
            def key(barcode):
                return hashlib.sha256(f'PDL-holdout-v1|{lib}|{barcode}'.encode()).digest()
            for index, barcode in enumerate(sorted(barcodes, key=key)):
                membership.append(dict(library=lib, cell_group=group, barcode=barcode, block=index % 4))
    assert len(membership) == 62490
    write_table(HERE / 'heldout_assignments.tsv.gz', membership)


def build(raw, features, cache, historical, check_raw):
    genes = read_rows(features)
    names = np.asarray([r['gene'] for r in genes])
    panel = np.asarray([r['ecm_component'] == 'other' for r in genes])
    assert len(set(names)) == len(names) == 10806 and sum(panel) == 101
    membership = read_rows(HERE / 'heldout_assignments.tsv.gz')
    records, arrays, depths, cells = [], [], [], []
    expected = None
    if historical:
        expected = np.load(historical, allow_pickle=False)
        lookup = {g:i for i,g in enumerate(expected['genes'])}
        rows_expected = [lookup[g] for g in names]
        assert np.array_equal(expected['classes'], CLASSES)
    published = {}
    if check_raw:
        published = {r['file']:r['sha256'] for r in read_rows(check_raw)}
    for lib in LIBS:
        paths = [next(raw.glob(f'*_{lib}_{f}')) for f in ('features.tsv.gz', 'barcodes.tsv.gz', 'matrix.mtx.gz')]
        if published:
            for path in paths:
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                assert digest == published[path.name]
        indices = {}
        for i,r in enumerate(read_rows(paths[0],header=False)):
            indices[r[1]] = i if r[1] not in indices else None
        take = [indices[g] for g in names]
        assert None not in take
        bars = [r[0] for r in read_rows(paths[1],header=False)]
        positions = {s:i for i,s in enumerate(bars)}
        assert len(bars) == len(positions)
        if cache:
            matrix = sparse.load_npz(cache / f'{lib}.npz').tocsr()
        else:
            with gzip.open(paths[2],'rb') as stream:
                matrix = io.mmread(stream).tocsr()
        rows = [r for r in membership if r['library'] == lib]
        ids = [positions[r['barcode']] for r in rows]
        assert len(set(ids)) == len(ids)
        columns = [CLASSES.index(r['cell_group'])*4 + int(r['block']) for r in rows]
        selector = sparse.csr_matrix((np.ones(len(ids), dtype=np.int64), (ids,columns)),
                                     shape=(matrix.shape[1], len(CLASSES)*4))
        summed = (matrix @ selector).toarray().reshape(matrix.shape[0],len(CLASSES),4)
        dat = summed[take].transpose(2,1,0).copy()
        depth = summed.sum(axis=0).T
        n = np.asarray(selector.sum(axis=0)).reshape(len(CLASSES),4).T
        assert np.all(n>=9) and np.all(n.max(axis=0)-n.min(axis=0)<=1)
        assert dat.min()>=0 and np.all(dat.sum(axis=2)<=depth)
        if expected is not None:
            np.testing.assert_array_equal(dat.sum(axis=0).T, expected[lib][rows_expected])
            np.testing.assert_array_equal(depth.sum(axis=0), expected[lib+'_total'])
        cells.append(n); depths.append(depth); arrays.append(dat)
        for block in range(4):
            for i,group in enumerate(CLASSES):
                records.append(dict(library=lib, block=block, cell_group=group,
                                    cells=int(n[block,i]), all_features_umi=int(depth[block,i])))
    if expected is not None:
        expected.close()
    assert sum(r['cells'] for r in records) == 62490
    write_table(HERE/'block_sizes.tsv', records)
    np.savez_compressed(HERE/'block_counts.npz', libraries=LIBS, groups=CLASSES, genes=names,
                        panel=panel, counts=np.stack(arrays).astype(np.int32),
                        total_umi=np.stack(depths), cells=np.stack(cells))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, default=WORK/'stage4/input')
    p.add_argument('--features', type=Path, default=WORK/'stage5/reference/features.tsv')
    p.add_argument('--cache', type=Path, help='Optional previously generated sparse matrices')
    p.add_argument('--historical', type=Path, help='Optional full-group counts for exact check')
    p.add_argument('--new-membership-from', type=Path, help='Regenerate saved assignments from old QC tables')
    p.add_argument('--check-raw',type=Path, help='Validate public file SHA-256 using stage4 manifest')
    args=p.parse_args()
    if args.new_membership_from:
        materialize_membership(args.raw,args.new_membership_from)
    build(args.raw,args.features,args.cache,args.historical,args.check_raw)


if __name__ == '__main__':
    main()
