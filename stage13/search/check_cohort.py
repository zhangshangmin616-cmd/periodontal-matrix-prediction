"""Verify the published Leuven sample design using public sources; write compact evidence."""
import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile

import pymupdf

HERE = Path(__file__).resolve().parent
THESIS = 'https://lirias.kuleuven.be/retrieve/0890d021-1dae-4af5-8def-ac15242d1730'
ARTICLE = 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12644353/fullTextXML'
SUPP = 'https://pmc-oa-opendata.s3.amazonaws.com/PMC12644353.1/'
ENA = 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1136245&result=read_run&fields=run_accession,sample_accession,sample_alias,experiment_alias,library_name,experiment_title'


def get(url):
    with urlopen(Request(url, headers={'User-Agent':'Mozilla/5.0'}), timeout=35) as f:
        return f.read()


def collapsed(text):
    return ' '.join(text.split())


def sample_tags(run):
    acc = run['sample_accession']
    e = ET.fromstring(get('https://www.ebi.ac.uk/ena/browser/api/xml/'+acc))
    sample = e.find('SAMPLE')
    assert sample.get('accession') == acc
    tags={a.findtext('TAG'): a.findtext('VALUE') for a in sample.findall('.//SAMPLE_ATTRIBUTE')}
    assert tags['isolate']==run['library_name'] and sample.get('alias')==run['sample_alias']
    return tags


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--document-cache',type=Path,help='Directory containing the unmodified public documents already downloaded')
    args=parser.parse_args()
    sources={'thesis':THESIS,'article':ARTICLE,'ENA':ENA}
    local={'thesis':'thesis.pdf','article':'secondary.xml','ENA':'source.xml'}
    for i in range(1,4):
        sources[f'supplement_{i}'] = SUPP+f'40510_2025_596_MOESM{i}_ESM.docx'
        local[f'supplement_{i}']=f'40510_2025_596_MOESM{i}_ESM.docx'
    with ThreadPoolExecutor(6) as pool:
        bodies={name:(args.document_cache/local[name]).read_bytes() for name in sources} if args.document_cache else dict(zip(sources,pool.map(get,sources.values())))
    with pymupdf.open(stream=bodies['thesis'],filetype='pdf') as pdf:
        tissue = collapsed(pdf[58].get_text())
        secondary = collapsed(pdf[114].get_text())
        assert len(pdf)==198 and '126 premolars' in tissue and 'combined in a single sample' in tissue
        assert '63 consented' in tissue and '46 and 41 patients' in secondary
    article = ET.fromstring(bodies['article'])
    abstract=collapsed(' '.join(article.find('.//abstract').itertext()))
    assert '69 DP and 63 PDL samples from 46 and 41 patients' in abstract
    text = collapsed(' '.join(article.find('.//body').itertext()))
    assert 'random effect for subject' in text
    for i in range(1,4):
        root=ET.fromstring(ZipFile(io.BytesIO(bodies[f'supplement_{i}'])).read('word/document.xml'))
        s=collapsed(' '.join(el.text or '' for el in root.iter('{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t')))
        if i==1:
            assert 'middle third of the root' in s
        else:
            assert 'Summary of findings' in s
    runs=list(csv.DictReader(io.StringIO(bodies['ENA'].decode()),delimiter='\t'))
    assert len(runs)==len({r['sample_accession'] for r in runs})==len({r['library_name'] for r in runs})==63
    with ThreadPoolExecutor(6) as pool:
        tags=list(pool.map(sample_tags,runs))
    keys=Counter(k for s in tags for k in s)
    donors=[k for k in keys if re.search(r'donor|participant|patient|subject|host_subject',k,re.I)]
    summary={'source_sha256':{name:{'url':sources[name],'sha256':hashlib.sha256(buf).hexdigest()}
                              for name,buf in bodies.items()},
             'dissertation_pdf_pages':{'bilateral_pooling_and_recruitment':59,'analysed_donors':115},
             'analysed_tissue_libraries':63,'analysed_pdl_donors_reported_in_secondary_paper':41,
             'registered_runs':len(runs),'unique_library_names':len({r['library_name'] for r in runs}),
             'sample_attributes_across_ENA_records':dict(sorted(keys.items())),
             'explicit_donor_fields':donors,'RNA_samples_by_treatment':dict(Counter(s['treatment'] for s in tags))}
    assert not donors and summary['RNA_samples_by_treatment']=={'Control':27,'7days':19,'28days':17}
    (HERE/'cohort_check.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8')
    print('Confirmed 63 sample records and published report of 41 PDL donors.')


if __name__=='__main__':
    main()
