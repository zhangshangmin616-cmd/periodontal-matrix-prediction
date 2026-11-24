"""Join a locally retrieved bibliography to its item-level eligibility decisions."""
import csv, html, json, re, sys
from pathlib import Path
from search import HERE, read, save


GROUPS = {
    "The force experiment used animals.": "11776918 18380969 22264327 25241549 25338791 25543591 27882108 31489407 32647179 33825762 37143931 37209709 38832299 40263606 41582843 41826577",
    "The force or other perturbation of the human cells was performed in culture.": "15256865 15562643 19555657 19602122 21120699 22332568 22389957 34185172 34692946 36055128",
    "Human cultures and rodent in vivo loading were studied in different parts of the paper.": "20578863 27161015 33892176 40763689",
    "The comparison concerned cell source, age, disease or tissue type and supplied no side-specific day-7 force contrast.": "22244633 22362338 25325708 28713908 31378921 35106969 36643271 38596952",
    "Review or hypothesis article.": "26599120 33923725 36220124 42874755",
    "Seven-day human physically separated pressure and tension PDL sampled; nine targeted RNA assays, one eligible matrix gene (TGFB1). Previously evaluated.": "17850423",
    "Post-loading extraction at 28 days. Measured ITGB1 and actin mRNA, without a day-7 arm.": "34041525",
    "Gingival tissue and gingival crevicular fluid were assayed.": "34693434",
    "Implant-site clinical case report with no pressure- and tension-specific day-7 PDL RNA.": "38742461",
    "Axial intrusion versus extrusion; days 0, 1, 3 and 10 and circumferential PDL collection.": "41706758",
    "Salivary biomarker measurements.": "41816270",
}


def main(path):
    meta={r["id"]:r for r in json.loads(path.read_text(encoding="utf8"))}
    hits=set(p for r in read(HERE/"pubmed_queries.tsv") for p in r["pmids"].split("|") if p)
    decisions={pm:reason for reason,ids in GROUPS.items() for pm in ids.split()}
    assert set(decisions)==set(meta)==hits
    result=[]
    for pm in sorted(hits,key=int):
        r=meta[pm]
        result.append(dict(pmid=pm,year=r.get("pubYear",""),title=html.unescape(re.sub("<[^>]*>","",r["title"])),
                           disposition=decisions[pm],source_url=f"https://pubmed.ncbi.nlm.nih.gov/{pm}/"))
    save(HERE/"screened.tsv",result)
    print("Screened",len(result),"distinct PubMed entries.")


if __name__=="__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE/"bibliography.json")
