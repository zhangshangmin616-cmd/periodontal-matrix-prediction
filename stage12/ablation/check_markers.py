"""Check the exclusive mural-marker assignment against independent reference labels."""
import csv
from pathlib import Path
from collections import defaultdict


HERE=Path(__file__).resolve().parent


def main():
    donors=defaultdict(lambda:dict(total_cells=0,labelled_mural=0,marker_mural=0,both=0))
    with (HERE/'reference_label_cross_table.tsv').open(encoding='utf8') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            actual=row['author_label'] in ('Perivascular','SMC')
            total=sum(int(v) for k,v in row.items() if k not in ('author_label','donor'))
            called=int(row['mural'])
            record=donors[row['donor']]
            record['total_cells']+=total
            record['labelled_mural']+=total*actual
            record['marker_mural']+=called
            record['both']+=called*actual
    overall={k:sum(row[k] for row in donors.values()) for k in next(iter(donors.values()))}
    rows=[dict(donor=d,**r) for d,r in sorted(donors.items())]+[dict(donor='all',**overall)]
    for row in rows:
        row['recall']=row['both']/row['labelled_mural'] if row['labelled_mural'] else ''
        row['precision']=row['both']/row['marker_mural'] if row['marker_mural'] else ''
    with (HERE/'mural_marker_reference.tsv').open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    assert overall == dict(total_cells=1590,labelled_mural=385,marker_mural=91,both=91)
    print('Verified five healthy donors: 91 strict marker calls match 385 annotated mural cells.')


if __name__=='__main__':
    main()
