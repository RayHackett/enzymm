import os
import json

def BB_score(Template):

    with open(Template, 'r') as f:
        Result = f.readlines()

    # We want to get a % of backbone RESIDUES in a template
    BB_dict = {}
    for line in Result:
        if line[:4] == 'ATOM':
            res = int(line[22:26])
            if res not in BB_dict:
                BB_dict[res] = 0 # False: is Sidechain
            if int(line[8:11]) >= 100:
                BB_dict[int(line[22:26])] = 1 # True: is Mainchain
    BB_score = sum(BB_dict.values())/len(BB_dict)
    
    return round(BB_score, 2) # just two digits after comma

with open('all_jess_templates_20230210.txt', 'r') as f:
    lines = f.readlines()
    
Template_list = []
for line in lines:
    Template_list.append(line.strip())

Template_BB_score_dict = {}
for Template in Template_list:
    Template_BB_score_dict[Template] = BB_score(Template)
    
with open('Template_BB_score_dict.json', 'w') as f:
    json.dump(Template_BB_score_dict, f)
