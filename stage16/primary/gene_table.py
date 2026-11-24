"""Make a one-row-per-gene table for the anatomically resolved RNA question."""
import csv
import gzip
from pathlib import Path
from search import HERE, read


def main():
    coverage = read(HERE / "coverage.tsv")
    with gzip.open(HERE.parents[1]/"stage10/state_multipliers.tsv.gz", "rt", encoding="utf8") as handle:
        effects = {r["gene"]: r for r in csv.DictReader(handle, delimiter="\t")}
    pooled = {r["gene"]:r for r in read(HERE.parent/"sides/forecasts.tsv") if r["assignment"] == "original"}
    assert {r["gene"] for r in coverage} == set(pooled)
    reference = effects["ACTB"]
    roles = {"biological target":"检测靶标", "reference only":"仅作内参",
             "no reported outcome in mapped assays":"未查到对应分侧RNA靶标结果"}
    rows = []
    for item in coverage:
        gene = item["gene"]
        effect, clinical = effects[gene], pooled[gene]
        rows.append({
            "基因":gene,
            "HGNC编号":item["hgnc_id"],
            "已核对的独立第7天人体分侧RNA证据":roles[item["role"]],
            "作为靶标的研究":item["measured_in"],
            "作为内参的研究":item["used_as_reference_in"],
            "压力侧相对未受力_组织实测方向":"上调" if item["pressure_vs_unloaded"] == "up" else "",
            "张力侧相对未受力_组织实测方向":"上调" if item["tension_vs_unloaded"] == "up" else "",
            "压力侧相对未受力_模型预测_log2相对ACTB":round(float(effect["PP7D"])-float(reference["PP7D"]), 6),
            "张力侧相对未受力_模型预测_log2相对ACTB":round(float(effect["TP7D"])-float(reference["TP7D"]), 6),
            "鲁汶未分侧组织_模型预测变化_百分位点":round(float(clinical["predicted_pooled"]), 6),
            "鲁汶未分侧组织_实测变化_百分位点":round(float(clinical["observed_unsplit_change"]), 6),
        })
    with (HERE/"gene_evidence.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
