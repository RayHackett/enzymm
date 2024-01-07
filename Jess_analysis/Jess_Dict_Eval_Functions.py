#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path

sys.path.append('/home/hackett/Software')
from common_functions import request_url
from common_functions import json_extract


with open(Path(Path(__file__).parents[0], 'template_res_dssp_dict.json'), 'r') as f:
    template_res_dssp_dict = json.load(f)

with open(Path(Path(__file__).parents[0], 'json_Enzymes_binding_sites.json'), 'r') as f:
    json_Enzymes_binding_sites = json.load(f)

with open(Path(Path(__file__).parents[0], 'json_NonEnzymes_binding_sites.json'), 'r') as f:
    json_NonEnzymes_binding_sites = json.load(f)

with open(Path(Path(__file__).parents[0], 'crystallization_hets.csv'), 'r') as f:
    lines = f.readlines()

cryst_artifacts = set()
for line in lines:
    if len(line.strip()) == 3:
        cryst_artifacts.add(line.strip())
cryst_artifacts.add(None) # some dont have the ligand annotated - we leave those out


# In[ ]:


def load_json(json_file):
    with open(json_file, 'r') as f:
        d = json.load(f)
    sort_dict(d)
    return d


# In[ ]:


def sort_dict(dictionary):
    myKeys = list(dictionary.keys())
    myKeys.sort(key=float)
    d = {i: dictionary[i] for i in myKeys}
    return d


# In[ ]:


def remove_result(dictionary):
    part_dict = dictionary
    full_dict = dictionary
    for key, key_values in full_dict.items():
        for sub_key, match_list in key_values.items():
            for i in range(len(match_list)):
                part_dict[key][sub_key][i].pop('Full_Result', '???')
                
    return part_dict


# In[ ]:


def filter_jess_dict(full_dict, value, operator, filter_value):
    part_dict = {}
    for key, key_values in full_dict.items():
        for sub_key, match_list in key_values.items():
            for match in match_list:
                expression = str(match[value]) + str(operator) + str(filter_value)
                if eval(expression):
                    if key not in part_dict:
                        part_dict[key] = {}
                    if sub_key not in part_dict[key]:
                        part_dict[key][sub_key] = []
                    part_dict[key][sub_key].append(match)

    return part_dict


# In[ ]:


def counter(dictionary):
    total_hits = 0
    keys = set()
    queries = set()
    templates = set()
    for key, key_values in dictionary.items():
        keys.add(key)
        for sub_key, match_list in key_values.items():
            queries.add(sub_key)
            for match in match_list:
                templates.add(match['Template'])
                total_hits += 1
    print('total hits, number of M-CSA entries, Number of Queries, Number of Templates')
    return total_hits, len(keys), len(queries), len(templates)


# In[ ]:


def count_resnum(dictionary):
    my_list = []
    for key, key_values in dictionary.items():
        for sub_key, match_list in key_values.items():
            for match in match_list:
                my_list.append(match['Temp_res_num'])

    res_dic = Counter(my_list)
    for key, value in res_dic.items():
        res_dic[key] = value
    
    template_sizes = [3, 4, 5, 6, 7, 8]
    for resnum in template_sizes:
        if resnum not in res_dic:
            res_dic[resnum] = 0
            
    res_dic = sort_dict(res_dic)
    
    return res_dic


# In[ ]:


def pie_chart(dictionary):
    
    fonts = {
    # use LaTeX fonts in the plot
    "font.family": 'serif',
    # Use 10pt font in plots, to match 11pt font in document
    "axes.labelsize": 11,
    "font.size": 11,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
    }

    plt.rcParams.update(fonts)

    # figsize=(width,heigth)
    cm = 1/2.54  # centimeters in inches
    textwidth = 14.69785*cm
    
    # Data to plot
    labels = []
    sizes = []

    for x, y in dictionary.items():
        labels.append(x)
        sizes.append(y)

    # Plot
    fig, axs = plt.subplots(figsize=(textwidth*0.5,5*cm), dpi=300)
    plt.pie(sizes, labels=labels)

    plt.axis('equal')
    plt.show()


# In[ ]:


def top_resnum_hits(dictionary):
    # get only the matches with the highest res numbers for each query
    highest_res_num = {}
    for key, key_values in d.items():
        for sub_key, match_list in key_values.items():
            resnum_list = []
            for match in match_list:
                resnum_list.append(int(match['Temp_res_num']))
            highest = np.argwhere(resnum_list == np.amax(resnum_list)).flatten().tolist()

            for i in highest:
                match = match_list[i]
                if key not in highest_res_num:
                    highest_res_num[key] = {}
                if sub_key not in highest_res_num[key]:
                    highest_res_num[key][sub_key] = []
                highest_res_num[key][sub_key].append(match)
    
    # there may be multiple clusters. these are redundant.
    # select the cluster which contains the hit with the lowest RMSD (all are same resnum here)
    highest_res_num_single_cluster = {}
    for key, key_values in highest_res_num.items():
        for sub_key, match_list in key_values.items():
            lowest_rmsd = 999
            best_match = ''
            for match in match_list: # find the lowest RMSD
                if match['RMSD'] <= lowest_rmsd:
                    lowest_rmsd = match['RMSD']
                    best_match = match
            cluster_number = best_match['Template'].split('.')[-4].split('_')[1]
            for match in match_list: # find the hits which share a cluster number with the hit with the lowest RMSD
                if match['Template'].split('.')[-4].split('_')[1] == cluster_number:
                    if key not in highest_res_num_single_cluster:
                        highest_res_num_single_cluster[key] = {}
                    if sub_key not in highest_res_num_single_cluster[key]:
                        highest_res_num_single_cluster[key][sub_key] = []
                    highest_res_num_single_cluster[key][sub_key].append(match)
                    
    return highest_res_num_single_cluster


# In[ ]:


def add_loc_check(jess_dict, json_data):
    for key, key_vals in jess_dict.items():
        for sub_key, match_list in key_vals.items():      
            for match in match_list:
                matched_resids = set()
                for line in match['Full_Result'][1:-1]:
                    matched_resids.add(int(line[22:26]))
                matched_resids = list(matched_resids)
                
                catalytic_resids = []
                try:
                    data = json_data[sub_key[:4]]
                    for site in data[sub_key[:4].lower()]:
                        for molecule in site['ligand_residues']:
                            if molecule['chem_comp_id'] in cryst_artifacts:
                                continue # this skips the for loop below at the same indentation as the if
                            for res in site['site_residues']: 
                                if res['chain_id'] == sub_key[4]: # author chain id
                                #if res['struct_asym_id'] == sub_key[4] # given chain id
                                    catalytic_resids.append(res['author_residue_number'])
                except:
                    pass

                match['catalytic_resids'] = catalytic_resids
                match['matched_resids'] = matched_resids
                
                annotated_resids = set(catalytic_resids) | set(match['Active_site_resids']) | set(match['Binding_site_resids'])
                
                match['correct_location'] = False
                for i in matched_resids:
                    if i in annotated_resids:
                        match['correct_location'] = True             

    return jess_dict


# In[ ]:


def remove_promiscuity(jess_dict, templates_to_remove):
    part_dict = {}
    for key, key_vals in jess_dict.items():
            for sub_key, match_list in key_vals.items():      
                for match in match_list:
                    #if side_chain_resids >= 3:
                    if match['Template'] not in templates_to_remove:
                        if key not in part_dict:
                            part_dict[key] = {}
                        if sub_key not in part_dict[key]:
                            part_dict[key][sub_key] = []
                        part_dict[key][sub_key].append(match)

    return part_dict


# In[ ]:


def frequency_of_hitnumbers(dictionary): # number of hits per query
    Query_dict = {}
    
    for key, key_values in dictionary.items():
        for sub_key, match_list in key_values.items():
            if sub_key not in Query_dict:
                Query_dict[sub_key] = 0
            Query_dict[sub_key] += len(match_list)
            
    mylist = list(Query_dict.values())
            
    return mylist


# In[ ]:


# are some templates too promisquious --> how many hits per template are there
def frequency_of_templatehits(dictionary):
    temp_dict = {} # dict with template names as keys and number of hits per template as value
    kies = set() # mcsa entries
    for key, key_values in dictionary.items():
        kies.add(key)
        for sub_key, match_list in key_values.items():
            for match in match_list:
                if match['Template'] not in temp_dict:
                    temp_dict[match['Template']] = 0
                temp_dict[match['Template']] += 1
            
    mylist = list(temp_dict.values())
            
    return mylist, temp_dict, kies


# In[ ]:


def frequency_of_hitnumbers_greater_3(dictionary):
    Query_dict = {}
    
    for key, key_values in dictionary.items():
        for sub_key, match_list in key_values.items():
            if sub_key not in Query_dict:
                Query_dict[sub_key] = 0
            for match in match_list:
                if match['Temp_res_num'] >3:
                    Query_dict[sub_key] += 1
            
    mylist = list(Query_dict.values())
            
    return mylist


# In[ ]:


def ec_match(match):
    ec_hit = False
    if match.get('Query_EC'): # is a list of EC numbers; evaluate each
        for ec in match['Query_EC']:
            if ec.split('.')[:3] == match['Template_EC'].split('.')[:3]:
                ec_hit = True
        return ec_hit
    else:
        return 'NA'
    

def loc_match(match):
    if match.get('Active_site_resids') or match.get('Binding_site_resids') or match.get('Catalytic_resids'):
        if match['correct_location'] == True:
            loc = True
        else:
            loc = False
        return loc
    else:
        return 'NA'

def cath_match(match):
    cath_hit = False
    if match.get('Query_CATH'): # is a list of EC numbers; evaluate each
        for cath in match['Query_CATH']:
            for t_cath in match['Template_CATH']:
                if cath.split('.')[:3] == t_cath.split('.')[:3]:
                    cath_hit = True
        return cath_hit
    else:
        return 'NA'
    

def interpro_match(match):
    if match['Query_InterPro'].get('domains') or match['Query_InterPro'].get('families') or match['Query_InterPro'].get('superfamilies'):
        if set(match['Template_InterPro']['domains']) & set(match['Query_InterPro']['domains']):
            hit = True
        elif set(match['Template_InterPro']['families']) & set(match['Query_InterPro']['families']):
            hit = True
        elif set(match['Template_InterPro']['superfamilies']) & set(match['Query_InterPro']['superfamilies']):
            hit = True
        else:
            hit = False
        return hit
    else:
        return 'NA'
    
    
def tester(match):
    max_score = 1
    score = 0
    if ec_match(match) != 'NA':
        max_score += 1
        if ec_match(match):
            score += 1
            
    if loc_match(match) != 'NA':
        max_score += 1
        if loc_match(match):
            score += 1
            
    if cath_match(match) != 'NA':
        max_score += 1
        if cath_match(match):
            score += 1
        
    if interpro_match(match) != 'NA':
        max_score += 1
        if interpro_match(match):
            score += 1
            
    if match['Preserved_resid_order']:
        score += 1
        
    return score/max_score

def calculate_similarity_in_rsas(Template, Full_Result):
    
    query_rsa_vals = []
    seen_resnums = set()
    for line in Full_Result[1:-1:3]:
        resname = line[17:20]
        chain = line[20:22].strip()
        resnum = int(line[22:26])
        if resnum not in seen_resnums:
            query_rsa_vals.append(float(line[55:60]))
        seen_resnums.add(resnum)

    if template_res_dssp_dict[Template]:
        template_rsa_vals = list(template_res_dssp_dict[Template].values())
    else:
        print('lacking rsas for', Template)
        return 'NA'

    rsa_differences = []
    if len(template_rsa_vals) == len(query_rsa_vals):
        for i in range(len(query_rsa_vals)):
                if query_rsa_vals[i] != 9.99 and template_rsa_vals[i] != 9.99:
                    rsa_differences.append(abs(query_rsa_vals[i]-template_rsa_vals[i]))
                    
    if len(rsa_differences) >= 3:
        return np.sqrt(np.mean(np.square(np.array(rsa_differences))))
    else:
        return float('NaN')


# In[ ]:


def create_cath_dict(d):
    cath_match_dict = {}
    full_match = 0
    level_3_match = 0
    level_2_match = 0
    level_1_match = 0
    no_match = 0
    no_data = 0
    total_matches = 0
    for key, key_vals in d.items():
        for sub_key, match_list in key_vals.items():
            for match in match_list:
                full_hit = False
                l3_hit = False
                l2_hit = False
                l1_hit = False
                no_hit = False
                total_matches += 1
                if match.get('Template_CATH'):
                    for Template_CATH in match['Template_CATH']:
                        if match.get('Query_CATH'): # is a list of EC numbers; evaluate each
                            for cath in match['Query_CATH']:
                                if cath == Template_CATH:
                                    full_hit = True
                                elif cath.split('.')[:3] == Template_CATH.split('.')[:3]:
                                    l3_hit = True
                                elif cath.split('.')[:2] == Template_CATH.split('.')[:2]:
                                    l2_hit = True
                                elif cath.split('.')[:1] == Template_CATH.split('.')[:1]:
                                    l1_hit = True
                                elif cath.split('.')[:1] != Template_CATH.split('.')[:1]:
                                    no_hit = True

                # for each of the Query_EC numbers only take the highest level of match
                if full_hit:
                    full_match += 1
                elif l3_hit:
                    level_3_match += 1
                elif l2_hit:
                    level_2_match += 1
                elif l1_hit:
                    level_1_match += 1
                elif no_hit:
                    no_match += 1
                else:
                    no_data += 1

    cath_match_dict['4'] = full_match
    cath_match_dict['3'] = level_3_match
    cath_match_dict['2'] = level_2_match
    cath_match_dict['1'] = level_1_match
    cath_match_dict['no'] = no_match
    cath_match_dict['NaN'] = no_data
    
    return cath_match_dict


# In[ ]:


def create_ec_dict(d):
    ec_dict = {}
    
    full_match = 0
    level_3_match = 0
    level_2_match = 0
    level_1_match = 0
    no_match = 0
    no_data = 0
    total_matches = 0
    unique_fullmatch_ec = set()
    for key, key_vals in d.items():
        for sub_key, match_list in key_vals.items():
            for match in match_list:
                full_hit = False
                l3_hit = False
                l2_hit = False
                l1_hit = False
                no_hit = False
                total_matches += 1
                Template_EC = match['Template_EC']
                if match.get('Query_EC'): # is a list of EC numbers; evaluate each
                    for ec in match['Query_EC']:
                        if pd.notna(ec):
                            if ec == Template_EC:
                                full_hit = True
                            elif ec.split('.')[:3] == Template_EC.split('.')[:3]:
                                l3_hit = True
                            elif ec.split('.')[:2] == Template_EC.split('.')[:2]:
                                l2_hit = True
                            elif ec.split('.')[:1] == Template_EC.split('.')[:1]:
                                l1_hit = True
                            elif ec.split('.')[:1] != Template_EC.split('.')[:1]:
                                no_hit = True

                # for each of the Query_EC numbers only take the highest level of match
                if full_hit:
                    full_match += 1
                elif l3_hit:
                    level_3_match += 1
                elif l2_hit:
                    level_2_match += 1            
                elif l1_hit:
                    level_1_match += 1                    
                elif no_hit:
                    no_match += 1
                else:
                    no_data += 1

    ec_dict['4'] = full_match
    ec_dict['3'] = level_3_match
    ec_dict['2'] = level_2_match
    ec_dict['1'] = level_1_match
    ec_dict['none'] = no_match
    ec_dict['NaN'] = no_data
    
    return ec_dict


# In[ ]:


def get_json_dict(jess_dict):
    json_dict = {}
    queries = set()
    for key, key_vals in jess_dict.items():
        for sub_key, match_list in key_vals.items():
            queries.add(sub_key)
            
    for query in queries:
        url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/binding_sites/' + query[:4]
        r = request_url(url, [200])
        if r:
            data = r.json()
            json_dict[query[:4]] = data

    return json_dict
            
#json_Enzymes_binding_sites = get_json_dict(Enzymes)
#json_NonEnzymes_binding_sites = get_json_dict(NonEnzymes)

#with open('json_Enzymes_binding_sites.json', 'w') as f:
#    json.dump(json_Enzymes_binding_sites, f)
#with open('json_NonEnzymes_binding_sites.json', 'w') as f:
#    json.dump(json_NonEnzymes_binding_sites, f)

