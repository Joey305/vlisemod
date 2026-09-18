import csv
import json
import re
from difflib import SequenceMatcher
from pathlib import Path


ROOT = Path('/Users/jxs794/Library/CloudStorage/OneDrive-UniversityofMiami/VLiSEMOD/Supplementary/Re-Submission-Package')
INPUT = Path('/tmp/citation_audit.json')
OUTPUT = ROOT / 'CITATION_SOURCE_CROSSWALK.csv'


def norm(value):
    return re.sub(r'[^a-z0-9]+', ' ', (value or '').lower()).strip()


def bib_numbers(report):
    result = {}
    for entry in report['bibliography_paragraphs']:
        match = re.match(r'\s*(\d+)\.?\s*(.*)', entry, re.S)
        if match:
            result[int(match.group(1))] = {'raw': match.group(2), 'norm': norm(match.group(2))}
    return result


def match_number(source, bibliography):
    title = source['title_norm']
    candidates = []
    doi = source.get('doi', '').lower()
    if doi:
        for number, entry_record in bibliography.items():
            if doi in entry_record['raw'].lower():
                return number
    for number, entry_record in bibliography.items():
        entry = entry_record['norm']
        entry_dois = re.findall(r'10\.\d{4,9}/[^\s)]+', entry_record['raw'].lower())
        if doi and entry_dois and doi not in entry_dois:
            continue
        score = SequenceMatcher(None, title, entry).ratio()
        if title and title in entry:
            score += 1
        candidates.append((score, number))
    score, number = max(candidates, default=(0, ''))
    return number if score >= 0.60 else ''


def action_for(doi, title, old_no, clean_no, markup_no):
    if doi == '10.1016/j.coph.2016.07.002':
        return 'REMOVE Southan from all citation groups and bibliography'
    if doi == '10.5281/zenodo.591637':
        return 'MERGE into version-specific RDKit source and REMOVE generic record'
    if doi == '10.5281/zenodo.19922430':
        return 'RETAIN as sole RDKit source and synchronize into marked manuscript'
    if doi in {'10.1021/acs.jcim.0c00589','10.1021/acs.jmedchem.8b01413','10.1016/j.bmcl.2024.129676','10.1038/s41467-018-08027-7'}:
        return 'RESTORE at opening Introduction second citation group'
    if doi in {'10.1021/acs.bioconjchem.9b00195','10.1111/j.1747-0285.2009.00921.x'}:
        return 'ADD at retrospective HIV-protease paragraph; replace static 45,46'
    if doi == '10.1101/2023.09.29.560163':
        return 'OMIT unless a current sentence requires citation'
    if doi == '10.1038/s41589-020-00689-z':
        return 'RETAIN one source; synchronize into marked and S5.7'
    return 'RETAIN; final number must be regenerated in Mendeley'


def status_for(doi, old_no, clean_no, markup_no):
    if doi == '10.1016/j.coph.2016.07.002': return 'REMOVE'
    if doi == '10.5281/zenodo.591637': return 'DUPLICATE/OBSOLETE'
    if doi in {'10.1021/acs.jcim.0c00589','10.1021/acs.jmedchem.8b01413','10.1016/j.bmcl.2024.129676','10.1038/s41467-018-08027-7','10.1021/acs.bioconjchem.9b00195','10.1111/j.1747-0285.2009.00921.x'}: return 'LOST'
    if doi == '10.1101/2023.09.29.560163': return 'REMOVE/UNUSED'
    if clean_no and markup_no: return 'PRESENT'
    if clean_no and not markup_no: return 'CLEAN ONLY'
    if markup_no and not clean_no: return 'MARKUP ONLY'
    return 'UNRESOLVED'


def main():
    reports = json.loads(INPUT.read_text())
    clean = reports['ACS_Infectious_Disease_2026_Schulz_etAL_Revised.docx']
    markup = reports['ACS_Infectious_Disease_2026_Schulz_etAL_MARKUP.docx']
    old = reports['ACS_Infectious_Disease_2026_Schulz_etAL_Reviewer_Pass11_FINAL-JMS.docx']
    docs = [(old, 'old'), (clean, 'clean'), (markup, 'markup')]
    bibs = {key: bib_numbers(report) for report, key in docs}
    sources = {}
    locations = {}
    for report, key in docs:
        for citation in report['citations']:
            for source in citation['items']:
                ident = source['identity']
                sources.setdefault(ident, source)
                if key == 'clean':
                    locations.setdefault(ident, []).append(f"clean P{citation['paragraph']}: {citation['context'][:90]}")
    for doi, title in [
        ('10.1021/acs.bioconjchem.9b00195','Acquired Removability of Aspartic Protease Inhibitors by Direct Biotinylation'),
        ('10.1111/j.1747-0285.2009.00921.x','How Much Binding Affinity Can be Gained by Filling a Cavity?'),
    ]:
        sources[doi] = {'doi': doi, 'title': title, 'title_norm': norm(title), 'identity': doi}
        locations[doi] = ['clean retrospective HIV-protease paragraph and S5.7']
    rows=[]
    for ident, source in sources.items():
        old_no=match_number(source,bibs['old'])
        clean_no=match_number(source,bibs['clean'])
        markup_no=match_number(source,bibs['markup'])
        doi=source['doi']
        rows.append({
            'source DOI': doi,
            'normalized title': source['title_norm'],
            'old Pass-11 number': old_no,
            'current clean number': clean_no,
            'current markup number': markup_no,
            'current manuscript citation locations': ' | '.join(locations.get(ident, [])),
            'current S5 citation locations': 'S5.6/S5.7 where applicable' if doi in {'10.1093/nar/28.1.235','10.1093/bioinformatics/btu789','10.1093/nar/gkad1004','10.1093/nar/gkv315','10.1016/j.jmb.2016.12.004','10.12688/f1000research.7931.1','10.1093/nar/gkae768','10.1093/nar/gkae1091','10.1021/acs.bioconjchem.9b00195','10.1111/j.1747-0285.2009.00921.x','10.1038/s41589-020-00689-z'} else '',
            'status': status_for(doi,old_no,clean_no,markup_no),
            'required action': action_for(doi,source['title'],old_no,clean_no,markup_no),
        })
    fields=list(rows[0])
    with OUTPUT.open('w', newline='') as handle:
        writer=csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(sorted(rows,key=lambda x:(x['old Pass-11 number']== '',int(x['old Pass-11 number'] or 999),x['source DOI'])))
    print(f'Wrote {len(rows)} source identities to {OUTPUT}')


if __name__ == '__main__':
    main()
