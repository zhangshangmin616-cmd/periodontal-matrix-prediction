"""Intersect published RNA targets with the previously fixed gene panel."""
from search import HERE, read, save


def main():
    assays = read(HERE / "assays.tsv")
    panel = read(HERE / "query_terms.tsv")
    rows = []
    for target in panel:
        gene = target["gene"]
        outcome = [a["study"] for a in assays if gene in a["biological_targets"].split("|")]
        reference = [a["study"] for a in assays if gene in a["reference_assayed"].split("|")]
        if outcome:
            role = "biological target"
        elif reference:
            role = "reference only"
        else:
            role = "no reported outcome in mapped assays"
        rows.append(dict(gene=gene, hgnc_id=target["hgnc_id"], role=role,
                         measured_in="|".join(outcome), used_as_reference_in="|".join(reference),
                         pressure_vs_unloaded="up" if gene == "TGFB1" else "",
                         tension_vs_unloaded="up" if gene == "TGFB1" else "",
                         numerical_pressure_minus_tension=""))
    assert len({r["gene"] for r in rows}) == len(rows) == 101
    assert [r["gene"] for r in rows if r["measured_in"]] == ["TGFB1"]
    assert [r["gene"] for r in rows if r["used_as_reference_in"]] == ["ACTB"]
    save(HERE / "coverage.tsv", rows)


if __name__ == "__main__":
    main()
