from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
import sys

HERE = Path(__file__).resolve().parent
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312900/"
FILES = {
    "geo_matrix.txt.gz": FTP+"matrix/GSE312900_series_matrix.txt.gz",
    "fpkm.xlsx": FTP+"suppl/GSE312900_fpkm.xlsx",
    "family.xml.tgz": FTP+"miniml/GSE312900_family.xml.tgz",
}

def get(entry):
    name,url = entry
    try:
        with urlopen(url,timeout=80) as f:
            data=f.read()
        if data[:2] not in [b"PK", bytes([31,139])]:
            return name,"Unexpected file format",data[:100]
        (HERE/name).write_bytes(data)
        return name,len(data)
    except OSError as e:
        return name,str(e)

if __name__=="__main__":
    entries = [(key,url) for key,url in FILES.items() if not sys.argv[1:] or key in sys.argv[1:]]
    with ThreadPoolExecutor(4) as executor:
        print(list(executor.map(get,entries)))
