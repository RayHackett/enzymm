import requests
import time
import json


mcsa_references = {}
for page in range(1,12):
    general = requests.get("https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page="+str(page))
    counter = 0
    while general.status_code != 200 and counter < 30:
        counter += 1
        time.sleep(1)
        general = requests.get("https://www.ebi.ac.uk/thornton-srv/m-csa/api/entries/?format=json&page="+str(page))
    if general.status_code == 200:
        general_page = general.json()
    else:
        print('page number' + page + 'does not exist:   ', general.status_code)

    for i in range(len(general_page['results'])):
        mcsa_references[general_page['results'][i]['mcsa_id']] = general_page['results'][i]['reference_uniprot_id']

def get_interpro(Uniprot_ID):
    interpro_domains = set()
    interpro_families = set()
    interpro_superfam = set()
    url = 'https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/' + Uniprot_ID
    r = requests.get(url)
    counter = 0
    while r.status_code not in [200, 204] and counter < 2:
        counter += 1
        time.sleep(1)
        r = requests.get(url)
        print(Uniprot_ID)
        print(r.status_code)
    if r.status_code == 200:
        data = r.json()

        for i in data['results']:
            if i['metadata']['type'] == 'domain':
                interpro_domains.add(i['metadata']['accession'])
            elif i['metadata']['type'] == 'family':
                interpro_families.add(i['metadata']['accession'])
            elif i['metadata']['type'] == 'homologous_superfamily':
                interpro_superfam.add(i['metadata']['accession'])

    return list(interpro_domains), list(interpro_families), list(interpro_superfam)

MCSA_interpro_dict = {}

for key, reference_seqs in mcsa_references.items():
    reference_seqs = reference_seqs.split(',')
    template_interpro_domains = set()
    template_interpro_fams = set()
    template_interpro_superfams = set()
    for seq in reference_seqs:
        seq = seq.replace(' ','')
        domains , fams, superfams = get_interpro(seq)
        template_interpro_domains.update(domains)
        template_interpro_fams.update(fams)
        template_interpro_superfams.update(superfams)

    MCSA_interpro_dict[key] = {'domains': list(template_interpro_domains), 'families': list(template_interpro_fams), 'superfamilies': list(template_interpro_superfams)}



with open('MCSA_interpro_dict.json', 'w') as f:
    json.dump(MCSA_interpro_dict, f)
