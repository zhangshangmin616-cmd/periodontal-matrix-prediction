"""Make an Excel-readable table of observed and predicted individual-gene rank changes."""
import csv
from pathlib import Path
from statistics import mean
HERE = Path(__file__).resolve().parent


def read(relative):
    with (HERE/relative).open(encoding="utf8",newline="") as f:
        return list(csv.DictReader(f,delimiter="\t"))


def main():
    wuhan = {(int(r["hours"]),r["gene"]):r for r in read("wuhan_counts/targets.tsv")
             if r["scale"]=="rank" and r["model"]=="tissue"}
    zhejiang = {r["gene"]:r for r in read("zhejiang/targets.tsv")
                if r["input"]=="NCBI_counts" and r["method"]=="tissue"}
    seoul = read("clinical_search/gene_changes.tsv")
    out = []
    for gene in sorted(zhejiang):
        row = {"基因":gene}
        for hour in [24,48]:
            for k,name in [("observed","观察"),("predicted","预测")]:
                row[f"武汉{hour}小时_{name}秩变化_百分位点"] = round(float(wuhan[hour,gene][k]),6)
        for k,name in [("observed","观察"),("predicted","预测")]:
            row[f"GSE106887_12小时_{name}秩变化_百分位点"] = round(float(zhejiang[gene][k]),6)
            vals=[float(r[k]) for r in seoul if r["gene"]==gene]
            assert len(vals)==10
            row[f"首尔自然咬合_10个配对组的平均{name}秩差_百分位点"] = round(mean(vals),6)
        out.append(row)
    assert len(out)==101 and len(wuhan)==202
    with (HERE/"gene_results.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(out[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)


if __name__ == "__main__":
    main()
