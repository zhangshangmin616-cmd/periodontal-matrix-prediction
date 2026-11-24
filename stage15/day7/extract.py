"""Re-extract the Chen human qPCR bar tops from two public figure sources."""
import argparse
import csv
from hashlib import sha256
import io
from pathlib import Path
import urllib.request
from zipfile import ZipFile

import numpy as np
import pymupdf
from PIL import Image


HERE = Path(__file__).resolve().parent
BASE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC10788340.1/"
SOURCES = {
    "PMC10788340.1.pdf": "e92b51431ccc18726f48b6e798c9840dea92565cf13f31f624172edda454ea48",
    "41368_2023_268_MOESM1_ESM.docx": "a95d076e3ca4c47941a1365783d3fc9ac03838f048a8f1545ba4c26e031c2c91",
}


def drawings(path):
    with pymupdf.open(path) as document:
        objects = document[5].get_drawings()
    for gene, x_left, x_right in [("TNFSF11", 200, 226), ("TNFRSF11B", 277, 304)]:
        bars = []
        for item in objects:
            rect = item["rect"]
            if (x_left < rect.x0 < rect.x1 < x_right and 247 < rect.y1 < 248
                    and 8 < rect.width < 9 and item["fill"] == (1., 1., 1.)
                    and item["items"][0][0] == "re"):
                bars.append(rect)
        control, loaded = sorted(bars, key=lambda rect: rect.x0)
        assert abs(control.y1 - loaded.y1) < .002
        yield [gene, "4d", "PDF point", round(control.y1, 3),
               round(control.y0, 3), round(loaded.y0, 3), .25]


def horizontal_centres(path):
    with ZipFile(path) as archive:
        source = archive.read("word/media/image3.tiff")
    assert sha256(source).hexdigest() == "8b68da089021a147b3b81476c16b614beaf7807113ab8b65a243744e64a69c96"
    img = np.asarray(Image.open(io.BytesIO(source)).convert("RGB"), dtype=int)
    assert img.shape == (3780, 2079, 3)
    panels = [("NLRP3", 554, 622, 618), ("CASP1", 797, 864, 611),
              ("GSDMD", 1054, 1119, 604), ("IL1B", 1304, 1368, 601)]
    for gene, control_x, loaded_x, axis_region in panels:
        dark = img[axis_region-4:axis_region+4, control_x-30:loaded_x+20].max(axis=2) < 70
        axis_rows = np.flatnonzero(dark.mean(axis=1) > .97) + axis_region-4
        assert len(axis_rows) in (2, 3) and np.all(np.diff(axis_rows) == 1)
        tops = []
        for x, colour in [(control_x, "blue"), (loaded_x, "magenta")]:
            red, green, blue = img[320:axis_region-4, x-16:x+16].transpose(2, 0, 1)
            match = ((green-red > 30) & (blue-red > 40) if colour == "blue"
                     else (red-green > 37) & (blue-green > 37))
            # A 31/32-pixel span captures the bar edge and excludes the narrower error cap.
            rows = np.where(match.sum(axis=1) >= 31)[0] + 320
            assert 2 <= len(rows) <= 5 and np.all(np.diff(rows) == 1)
            tops.append(float(np.mean(rows)))
        yield [gene, "S3a", "image pixel", float(axis_rows.mean()), *tops, 2.]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=HERE / "coordinates.tsv")
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, digest in SOURCES.items():
        path = args.cache / name
        if not path.is_file():
            urllib.request.urlretrieve(BASE + name, path)
        assert sha256(path.read_bytes()).hexdigest() == digest
        paths.append(path)
    rows = [*horizontal_centres(paths[1]), *drawings(paths[0])]
    with args.output.open("w", encoding="utf8", newline="") as stream:
        out = csv.writer(stream, delimiter="\t", lineterminator="\n")
        out.writerow(["gene", "figure", "units", "baseline_y", "control_top", "loaded_top", "tolerance"])
        out.writerows(rows)


if __name__ == "__main__":
    main()
