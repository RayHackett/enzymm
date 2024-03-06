import requests
import os
import time
import json


def main():
    ec_dict = {}
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
            print('page number' + str(page) + 'does not exist:   ', general.status_code)
 
        for i in range(len(general_page['results'])):
            mcsa_id = int(general_page['results'][i]['mcsa_id'])
            if mcsa_id not in ec_dict:
                ec_dict[mcsa_id] = ''
            if len(general_page['results'][i]['all_ecs']) > 1:
                print(mcsa_id)
            for ec_num in general_page['results'][i]['all_ecs']:
                ec_dict[mcsa_id] = ec_num
 
    with open('MCSA_EC_mapping.json', 'w') as outfile:
        json.dump(ec_dict, outfile)

    return ec_dict

if __name__ == "__main__":
    main()
