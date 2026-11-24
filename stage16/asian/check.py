"""Verify the saved repository search results without network access."""
import csv
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def rows(filename):
    with (HERE / filename).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main():
    records = rows("ddbj_records.tsv")
    queries = rows("ddbj_queries.tsv")
    controls = rows("ddbj_controls.tsv")
    assert len(records) == len({(r["type"], r["identifier"]) for r in records})
    assert Counter(r["type"] for r in records) == {"sra-study":149, "bioproject":217}
    assert len(queries) == 12
    memberships = Counter((r["type"], term)
                          for r in records for term in r["matched_queries"].split(" | "))
    for q in queries:
        assert memberships[(q["type"], q["term"])] == int(q["total"]), q
    assert {r["type"] for r in controls if int(r["total"]) > 0} == {"gea","jga-study"}
    ids = {r["identifier"] for r in records}
    assert {"PRJNA1450544", "PRJNA1136245", "PRJDB19639"} <= ids
    sources = rows("sources.tsv")
    assert len(sources) == len({r["id"] for r in sources})
    assert all(r["primary_url"].startswith("https://") for r in sources)
    print("Verified", len(records), "repository records,", len(queries),
          "searches,", len(controls), "controls and", len(sources), "primary-source decisions.")


if __name__ == "__main__":
    main()
