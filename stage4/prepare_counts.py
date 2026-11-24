"""Cache sparse counts and compute per-library quality flags."""
import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io, sparse

ROOT = Path(__file__).resolve().parent
INPUT, CACHE, RESULTS = (ROOT / key for key in ("input", "cache", "results"))
SAMPLES = ["CP", "PP3D", "PP7D", "TP3D", "TP7D"]
PANELS = {
    "fibroblast_like": ["MDK","SFRP4","PDGFRA","POSTN"],
    "endothelial": ["PECAM1","CDH5","CLDN5","PLVAP"],
    "mural": ["RGS5","NOTCH3","MYH11","KCNJ8"],
    "immune": ["PTPRC","CORO1A","LAPTM5","CD53"],
    "glial": ["SOX10","PLP1","S100B","GFRA3"],
    "epithelial": ["KRT5","KRT14","KRT19","EPCAM"],
    "erythroid": ["HBB","HBA1","HBA2","ALAS2"],
}

def source_files(sample):
    return [next(INPUT.glob(f"*_{sample}_{filename}")) for filename in
            ["features.tsv.gz", "barcodes.tsv.gz", "matrix.mtx.gz"]]

def feature_table(path):
    return pd.read_csv(path, sep="\t", header=None, names=["id","symbol","type"],
                       dtype=str, keep_default_na=False)

def bounds(values, log):
    y = np.log1p(values) if log else values
    mid = np.median(y)
    mad = 1.4826 * np.median(np.abs(y - mid))
    lo, hi = mid - 3 * mad, mid + (4.5 if log else 3) * mad
    if log:
        return max(1, int(np.ceil(np.expm1(lo)))), int(np.floor(np.expm1(hi)))
    return lo, max(5.0, hi)

def totals(mx, mask=None):
    return np.asarray((mx if mask is None else mx[mask]).sum(axis=0)).ravel()

def run(sample):
    CACHE.mkdir(exist_ok=True)
    RESULTS.mkdir(exist_ok=True)
    feature_path, barcode_path, matrix_path = source_files(sample)
    f = feature_table(feature_path)
    barcodes = pd.read_csv(barcode_path, header=None, dtype=str)[0]
    assert f.id.is_unique and barcodes.is_unique
    with gzip.open(matrix_path,"rb") as stream:
        mx = io.mmread(stream).astype(np.int32).tocsc()
    assert mx.shape == (len(f),len(barcodes)) and mx.data.min() > 0
    mx.sum_duplicates()
    mx.sort_indices()
    sparse.save_npz(CACHE/f"{sample}.npz",mx,compressed=False)
    umi = totals(mx)
    genes = np.diff(mx.indptr)
    mt = 100 * totals(mx,f.symbol.str.startswith("MT-").to_numpy()) / np.maximum(umi,1)
    hb = 100 * totals(mx,f.symbol.isin(["HBA1","HBA2","HBB","HBD","HBG1","HBG2"]).to_numpy()) / np.maximum(umi,1)
    gl, gh = bounds(genes, True)
    ul, uh = bounds(umi, True)
    _, mh = bounds(mt, False)
    pass_qc = (genes >= gl) & (genes <= gh) & (umi >= ul) & (umi <= uh) & (mt <= mh)
    qc = pd.DataFrame(dict(barcode=barcodes, umi=umi, genes=genes, mt_pct=mt, hb_pct=hb, pass_qc=pass_qc))
    qc.to_csv(CACHE/f"{sample}_qc.tsv.gz",sep="\t",index=False)
    qc.loc[pass_qc,["barcode"]].to_csv(CACHE/f"{sample}_dbl_barcodes.tsv",index=False,header=False)
    f.id.to_csv(CACHE/f"{sample}_dbl_genes.tsv",index=False,header=False)
    selected = mx[:,pass_qc]
    for name, values in (("i",selected.indices),("p",selected.indptr),("x",selected.data)):
        values.astype("<i4",copy=False).tofile(CACHE/f"{sample}_dbl_{name}.bin")
    # Compact binary CSC stores the same counts as a Matrix Market matrix without text parsing.
    rec = dict(library=sample,n_input=len(barcodes),n_qc=int(pass_qc.sum()),n_removed=int((~pass_qc).sum()),
               genes_low=gl,genes_high=gh,umi_low=ul,umi_high=uh,mt_high=mh,
               low_genes=int((genes<gl).sum()),high_genes=int((genes>gh).sum()),low_umi=int((umi<ul).sum()),
               high_umi=int((umi>uh).sum()),high_mt=int((mt>mh).sum()),median_genes=np.median(genes),
               median_umi=np.median(umi),median_mt=np.median(mt),median_hb=np.median(hb))
    pd.DataFrame([rec]).to_csv(CACHE/f"{sample}_thresholds.tsv",sep="\t",index=False)
    print(sample,len(barcodes),pass_qc.sum(),gl,gh,ul,uh,round(mh,2),flush=True)

if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("samples",nargs="*",default=SAMPLES)
    for sample in parser.parse_args().samples:
        run(sample)
    available=[CACHE/f"{key}_thresholds.tsv" for key in SAMPLES]
    if all(x.exists() for x in available):
        pd.concat([pd.read_csv(x,sep="\t") for x in available]).to_csv(RESULTS/"qc_thresholds.tsv",sep="\t",index=False)
