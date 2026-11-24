"""Rebuild the figure and supplementary-table data from the archived analysis inputs."""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def run(rel):
    subprocess.run([sys.executable, str(ROOT/rel)], cwd=ROOT, check=True)

def main():
    checks = json.loads((ROOT/"sha256.json").read_text(encoding="utf8"))
    for path, digest in checks.items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != digest:
            raise ValueError("Input checksum differs: "+path)
    if "--from-counts" in sys.argv:
        for task in ("stage17/wuhan_counts/score.py", "stage17/zhejiang/score.py",
                     "stage18/transfer/score.py","stage18/deviations/score.py",
                     "stage18/bounds/rank_limits.py","stage18/check/verify.py"):
            run(task)
    run("stage19/figures/assemble_data.py")
    run("stage19/supp/make_source_data.py")
    for path in sorted((ROOT/"stage19").rglob("*.csv")):
        key=path.relative_to(ROOT).as_posix()
        original=checks[key]
        if hashlib.sha256(path.read_bytes()).hexdigest()!=original:
            raise ValueError("Changed output: "+key)
    with (ROOT/"stage19/figures/Figure_1_source_data.csv").open(encoding="utf8",newline="") as f:
        rows=list(csv.DictReader(f))
    assert len(rows)==101
    closer=sum(abs(float(s["atlas"])-float(s["observed"])) < abs(float(s["observed"])) for s in rows)
    assert closer==41
    if "--figures" in sys.argv:
        run("stage19/figures/plot_figures.py")
    print("Rebuilt data for the three figures and seven supplementary tables; verified the seven-study design table.")

if __name__=="__main__":
    main()
