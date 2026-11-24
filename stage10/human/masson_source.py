"""Read the author's anatomical group medians and labelled force-region means."""
import csv
import hashlib
import math
from pathlib import Path
import re
from urllib.request import urlopen
from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
URL = 'https://repositorio.pgsscogna.com.br/bitstream/123456789/2736/1/282d74066b1278126458adace089b49c.pdf'
EXPECTED = 'd4c9701073be82719eb7043af8ddac7e67b7dbfb26c6ae4867f99ded7531133b'


def save(name, keys, records):
    with (HERE/name).open('w',encoding='utf-8',newline='') as f:
        writer = csv.writer(f,delimiter='\t',lineterminator='\n')
        writer.writerow(keys)
        writer.writerows(records)


def main():
    pdf = HERE/'masson.pdf'
    if not pdf.exists():
        with urlopen(URL, timeout=30) as res:
            pdf.write_bytes(res.read())
    assert hashlib.sha256(pdf.read_bytes()).hexdigest() == EXPECTED
    pages = PdfReader(pdf).pages
    text = pages[30].extract_text()
    # Original printed page 29. Medians followed by first and third quartiles.
    pattern = re.compile(r'(?m)^(VC|VM|VA|LA|LM|LC)\s+(.+)$')
    timedays = [0,3,7,15,21]
    records = []
    for i, (site,row) in enumerate(pattern.findall(text)):
        values = re.findall(r'(\d+,\d+)\s*\(\s*(\d+,\d+)\s*-\s*(\d+,\d+)\s*\)',row)
        assert len(values)==5, (site,row)
        for day,stats in zip(timedays,values,strict=True):
            median,q1,q3 = [float(x.replace(',','.')) for x in stats]
            assert q1 <= median <= q3
            records.append(['IFNG' if i<6 else 'IL10',site,day,median,q1,q3])
    assert len(records)==60
    save('masson_anatomical_medians.tsv', ['gene','anatomical_site','day','median','q1','q3'],records)
    # Page 30, Table 3, values transcribed directly from PDF. SD is descriptive.
    regions = [('IFNG',4.08,1.06,7.24,.95,5.50,1.21),
               ('IL10',3.92,1.25,6.37,1.44,4.99,1.06)]
    text = re.sub(r'\s+','',pages[31].extract_text()).replace(',','.')
    for _,t,td,p,pd,m,md in regions:
        for mean,sd in [(t,td),(p,pd),(m,md)]:
            assert f'{mean:.2f}±{sd:.2f}' in text
    save('masson_zone_day7.tsv',
         ['gene','day','reference','tension_mean','tension_sd','pressure_mean','pressure_sd',
          'intermediate_mean','intermediate_sd','log2_C_T','log2_C_T_print_rounding_low','log2_C_T_print_rounding_high'],
         [[g,7,'ACTB',t,td,p,pd,m,md,math.log2(p/t),math.log2((p-.005)/(t+.005)),
           math.log2((p+.005)/(t-.005))] for g,t,td,p,pd,m,md in regions])


if __name__=='__main__':
    main()
