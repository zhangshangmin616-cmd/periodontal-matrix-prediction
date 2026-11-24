"""Record the narrowly targeted public registry searches, independently of model errors."""
import csv
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

root = Path(__file__).resolve().parent


def fetch(url):
    with urlopen(Request(url, headers={'User-Agent': 'Mozilla/5.0'}), timeout=30) as response:
        return response.read()


def query(kind, words):
    if kind == 'ENA':
        base = 'https://www.ebi.ac.uk/ena/portal/api/search'
        params = dict(result='read_run', query=f'tax_eq(9606) AND library_strategy="RNA-Seq" AND study_title="{words}"',
                      fields='study_accession,run_accession,study_title', format='tsv', limit=10000)
    elif kind == 'ClinicalTrials.gov':
        base = 'https://clinicaltrials.gov/api/v2/studies'
        params = {'query.term': words, 'pageSize': 100, 'countTotal': 'true',
                  'fields': 'NCTId,BriefTitle,HasResults,OverallStatus'}
    elif kind == 'GEO':
        base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        params = dict(db='gds', term=words, retmax=250, retmode='json')
    else:
        raise ValueError(kind)
    url = base+'?'+urlencode(params)
    content=fetch(url)
    if kind=='ENA':
        result=list(csv.DictReader(io.StringIO(content.decode()),delimiter='\t'))
        count=len(result)
        if count >= 10000:
            raise ValueError("increase limit before interpreting the ENA output")
        groups={}
        for row in result:
            group=groups.setdefault(row['study_accession'],{'title':row['study_title'],'runs':0})
            group['runs']+=1
        result=groups
    elif kind=='GEO':
        raw=json.loads(content)['esearchresult']
        result=raw['idlist']
        count=int(raw['count'])
    else:
        raw=json.loads(content)
        count=raw['totalCount']
        result=[{'nct': x['protocolSection']['identificationModule']['nctId'],
                 'title': x['protocolSection']['identificationModule']['briefTitle'],
                 'results': x.get('hasResults',False)} for x in raw.get('studies',[])]
    return {'registry':kind, 'query':words, 'url':url, 'total':count, 'returned':result}


queries=[
    ('ENA','*orthodontic*'),('ENA','*periodontal*'),
    ('ClinicalTrials.gov','"periodontal ligament" AND (orthodontic OR "tooth movement") AND (transcriptome OR sequencing OR expression)'),
    ('GEO','("orthodontic" OR "tooth movement") AND "Homo sapiens"[Organism] AND gse[ETYP]'),
]
if __name__=='__main__':
    with ThreadPoolExecutor(4) as pool:
        calls=[(pair,pool.submit(query,*pair)) for pair in queries]
        results=[]
        for (kind,words),task in calls:
            try:
                record=task.result();print(kind,record['total'])
            except Exception as e:
                record={'registry':kind,'query':words,'error':str(e)};print(record)
            results.append(record)
        (root/'registry_hits.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf8')
        geo = next((record for record in results if record['registry']=='GEO' and 'returned' in record), None)
        if geo:
            url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?'+urlencode(dict(
                db='gds',id=','.join(geo['returned']),retmode='json'))
            meta=json.loads(fetch(url))['result']
            fields=['accession','title','summary','n_samples','gdsType','taxon']
            records=[{field:meta[key].get(field) for field in fields} for key in geo['returned']]
            (root/'geo_candidates.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf8')
