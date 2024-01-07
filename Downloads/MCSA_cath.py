import requests
import json
import time

mcsa_cath = {}
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
        cath_ids = set()
        for resid in general_page['results'][i]['residues']:
            for chain in resid['residue_chains']:
                cath_ids.add(chain['domain_cath_id'])

        if '' in cath_ids:
            cath_ids.remove('')
        cath_ids = list(cath_ids)
        mcsa_cath[general_page['results'][i]['mcsa_id']] = cath_ids
        
with open('/homes/hackett/Downloads/MCSA_CATH_mapping.json', 'w') as f:
    json.dump(mcsa_cath, f)
