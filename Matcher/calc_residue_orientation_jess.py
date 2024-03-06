import argparse
import numpy as np
import itertools
from scipy.spatial import distance

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
            
def get_vec_template(residue):
    # dictionary in which the vectors from start to finish are defined for each residue type
    vector_atom_type_dict = {'GLY': ('C','O'),
                         'ALA': ('CA','CB'),
                         'VAL': ('CA','CB'),
                         'LEU': ('CA','CB'),
                         'ILE': ('CA','CB'),
                         'MET': ('CG','SD'),
                         'PHE': ('CZ','mid'),
                         'TYR': ('CZ','OH'),
                         'TRP': ('CZ2','NE1'),
                         'CYS': ('CB','SG'),
                         'PRO': ('C','O'),
                         'SER': ('CB','OG'),
                         'THR': ('CB','OG1'),
                         'ASN': ('CG','OD1'),
                         'GLN': ('CD','OE1'),
                         'LYS': ('CE','NZ'),
                         'ARG': ('CZ','mid'),
                         'HIS': ('CG','ND1'),
                         'ASP': ('CG','mid'),
                         'GLU': ('CD','mid'),
                         'PTM': ('CA','CB'),
                         'ANY': ('C','O')}
    
    test_atom = list(residue)[0]
    res_type = test_atom[17:20].strip()
    vectup = vector_atom_type_dict[res_type]
    
    if vectup[1] != 'mid': # from first atom to second atom
        for index, atom in enumerate(residue):
            atom_name = atom[12:16].strip()
            if atom_name == vectup[0]:
                first_atom_index = index
                first_x = float(atom[30:38])
                first_y = float(atom[38:46])
                first_z = float(atom[46:54])
            elif atom_name == vectup[1]:
                second_atom_index = index
                second_x = float(atom[30:38])
                second_y = float(atom[38:46])
                second_z = float(atom[46:54])
                
        vec = [first_x-second_x, first_y-second_y, first_z-second_z]   
        return vec, [first_atom_index, second_atom_index], res_type
            
    else: # from the middle atom to the midpoint
        first_atom = False
        for index, atom in enumerate(residue):
            atom_name = atom[12:16].strip()
            if atom_name == vectup[0]:
                middle_atom_index = index
                middle_x = float(atom[30:38])
                middle_y = float(atom[38:46])
                middle_z = float(atom[46:54])
            elif first_atom == False:
                first_atom = True
                first_atom_index = index
                first_atom_name = atom_name
                side1_x = float(atom[30:38])
                side1_y = float(atom[38:46])
                side1_z = float(atom[46:54])
            elif first_atom == True:
                second_atom_index = index
                second_atom_name = atom_name
                side2_x = float(atom[30:38])
                side2_y = float(atom[38:46])
                side2_z = float(atom[46:54])

        # calculate midpoint between the two side atoms        
        midpoint = [(side1_x+side2_x)/2, (side1_y+side2_y)/2, (side1_z+side2_z)/2]
        # calculate the vector between middle atom and midpoint
        vec = [middle_x-midpoint[0], middle_y-midpoint[1], middle_z-midpoint[2]]

        return vec, [middle_atom_index, first_atom_index, second_atom_index], res_type

def get_query_vec(residue, vec_info):
    # we use Template_residue to know which atom was used as the middle atom for calculating the template mid-vector of that residue
    if len(vec_info) == 3:
        for index, atom in enumerate(residue):
            if index == vec_info[0]: # middle_atom:
                middle_x = float(atom[30:38])
                middle_y = float(atom[38:46])
                middle_z = float(atom[46:54])
            elif index == vec_info[1]: # side_atom 1
                side1_x = float(atom[30:38])
                side1_y = float(atom[38:46])
                side1_z = float(atom[46:54])
            elif index == vec_info[2]: # side_atom 2
                side2_x = float(atom[30:38])
                side2_y = float(atom[38:46])
                side2_z = float(atom[46:54])

        # calculate midpoint between the two side atoms        
        midpoint = [(side1_x+side2_x)/2, (side1_y+side2_y)/2, (side1_z+side2_z)/2]
        # calculate the vector between middle atom and midpoint
        vec = [middle_x-midpoint[0], middle_y-midpoint[1], middle_z-midpoint[2]]

    elif len(vec_info) == 2: # uses just two atoms
        for index, atom in enumerate(residue):
            if index == vec_info[0]: # first atom
                first_x = float(atom[30:38])
                first_y = float(atom[38:46])
                first_z = float(atom[46:54])
            elif index == vec_info[1]: # second atom
                second_x = float(atom[30:38])
                second_y = float(atom[38:46])
                second_z = float(atom[46:54])
                    
        vec = [first_x-second_x, first_y-second_y, first_z-second_z]
        
    return vec
        
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2): # from https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def calc_angle_rms(angle_list):
    """ Calculates the RMS angles"""
    return np.sqrt(np.mean(np.square(angle_list)))

def calc_angle_mean(angle_list):
    """simple mean of angles"""
    return np.mean(angle_list)
    
def calc_angle_rmsd(angle_list1, angle_list2):
    """ Calculates the RMSD of two lists of angles"""
    return np.sqrt(np.mean(np.square(np.array(angle_list1)-np.array(angle_list2))))

def open_jess_pdb(pdb_file):
    with open(pdb_file, 'r') as f:
        Result = f.readlines()
        
    template_atom_lines = [] # list of atom lines
    for line in Result:
        if line[:4] == 'ATOM': # get only the atom lines
            template_atom_lines.append(line)
    
    return chunks(template_atom_lines, 3) # generator that yields 3 lines equal to one residue at a time
        
def calculate_all_vectors(Template, Match):
    
    temp_residue_generator = open_jess_pdb(Template)
    template_vec_list = []
    vec_info_list = []
    for residue in temp_residue_generator: # residue order in the template and the match are the same
        # get_vec_template returns the format: vec, atom_index_list, res_type
        result = get_vec_template(residue)
        vec_info_list.append(result)
        template_vec_list.append(result[0]) # append the vector
        
    match_residue_generator = open_jess_pdb(Match)
    match_vec_list = []
    for i, residue in enumerate(match_residue_generator): # now iterate over matched residues
        vec_info = vec_info_list[i][1] # incoorporating information about atom_index_list from template
        match_vec = get_query_vec(residue, vec_info)
        match_vec_list.append(match_vec)
            
    return template_vec_list, match_vec_list

def angle_rms(temp_vec_list, match_vec_list):
    # now calculate the angle between the vector of the template and the query per residue
    angle_list = []
    for i in range(len(temp_vec_list)):
        angle_list.append(angle_between(match_vec_list[i], temp_vec_list[i]))

    # calculate the RMS of the angle list. each angle in angle_list is for one residue
    angle_rms_val = calc_angle_rms(angle_list)
    return angle_rms_val

def angle_mean(temp_vec_list, match_vec_list):
    # now calculate the angle between the vector of the template and the query per residue
    angle_list = []
    for i in range(len(temp_vec_list)):
        angle_list.append(angle_between(match_vec_list[i], temp_vec_list[i]))

    # calculate the RMS of the angle list. each angle in angle_list is for one residue
    angle_mean_val = calc_angle_mean(angle_list)
    return angle_mean_val

def calc_all_angle_permutations(vec_list):
    """ Calculate all pairwise permutations of angles from a list of vectors """
    # generate list of possible pairwise permutations of angles
    pair_order_list = itertools.combinations(vec_list,2)
    # now calculate the angles between the possible permutations
    angle_list = []
    for pair in pair_order_list:
        angle_list.append(angle_between(pair[0], pair[1]))
    return angle_list
        
def angle_difference_test(temp_vec_list, match_vec_list):
    # calculate the pairwise angles
    temp_angle_list = calc_all_angle_permutations(temp_vec_list)
    match_angle_list = calc_all_angle_permutations(match_vec_list)

    # calculate the rmsd between the two lists of angles
    angle_rmsd = calc_angle_rmsd(temp_angle_list, match_angle_list)
    return angle_rmsd

def get_residue_positions(pdb_file):
    residue_generator = open_jess_pdb(pdb_file)
    residue_positions = []
    for residue in residue_generator:
        res_pos = np.empty(shape=(3, 3)) # initialize an empty numpy array, rows are atoms, columns are x, y, z
        for index, atom in enumerate(residue): # iterate over atoms in residue
            atom_x = float(atom[30:38])
            atom_y = float(atom[38:46])
            atom_z = float(atom[46:54])
            res_pos[index] = np.array([atom_x, atom_y, atom_z]) # fill the array at a given row with coordinates of an atom
        residue_positions.append(res_pos)
        
    return residue_positions

# calculate all pairwise distances
def calculate_minimum_pairwise_distance(residue_positions):
    # generate list of possible pairwise residue distances
    residue_pair_list = itertools.combinations(residue_positions,2)

    pairwise_dists = []
    for pair in residue_pair_list:
        pairwise_dists.append(calc_min_dist_between_residues(pair[0], pair[1]))
        
    return pairwise_dists

def calc_min_dist_between_residues(res1, res2):
    """ res1 and res2 are each numpy arrays with XYZ coordinates of 3 atoms each - basically arrays of size 3,3 """
    if np.shape(res1) == (3,3) and np.shape(res2) == (3,3): # just a check if we have xzy coordinates for all 3 atoms
        # distance.cdist computes euclidian distances between each pair of the two residues
        return distance.cdist(res1,res2, metric='euclidean').min() # minimum of the np array
    else:
        print('Wrong array shape for residues')
        
def calculate_weighted_rmsd(angle_list1, angle_list2, weights):
    # weights should be the distances from the template!
    """ Rescale the distances into scaled weights - smaller distance is better"""
    max_val = 10
    min_val = 0.5
    # basically like an inverted min max scaling between 0 and 1 but instead between min_val and max_val
    # inverted since small distances are better and should be weighted higher
    scaled_weights = [min_val + (max(weights)-i)*(max_val-min_val)/(max(weights)-min(weights)) for i in weights]
    
    """Sum over (squared difference of angles)*weight divided by sum over weights"""
    return np.sqrt( np.sum(np.multiply(scaled_weights, np.square(np.array(angle_list1)-np.array(angle_list2)))) / np.sum(scaled_weights) )

def weighted_angle_difference_test(temp_vec_list, match_vec_list, template_residue_positions):
    # calculate the pairwise angles
    temp_angle_list = calc_all_angle_permutations(temp_vec_list)
    match_angle_list = calc_all_angle_permutations(match_vec_list)
    
    # calculate the minimum pairwise distances between residues
    template_distances = calculate_minimum_pairwise_distance(template_residue_positions)

    # calculate the rmsd between the two lists of angles
    weighted_angle_rmsd = calculate_weighted_rmsd(temp_angle_list, match_angle_list, template_distances)
    return weighted_angle_rmsd
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--template', help='Jess Template in .pdb format')
    parser.add_argument('-m', '--match', help='Jess Match with that Template in .pdb format')
    parser.add_argument('-avg', help='Simple arithmetic mean of Angle between corresponding Template and Query Residues', action='store_true', default=False)
    parser.add_argument('-rms', help='Root Mean Square of Angle between corresponding Template and Query Residues', action='store_true', default=False)
    parser.add_argument('-adt', help='Angle Difference Test for RMSD of all pairwise Angles in Template and Query', action='store_true', default=False)
    parser.add_argument('-wadt', help='Weighted Angle Difference Test by distances in the Template for RMSD of all pairwise Angles in Template and Query', action='store_true', default=False)
    args = parser.parse_args()
    
    Template = args.template
    Match = args.match
    
    # calculate the vectors of the template and of the match
    temp_vec_list, match_vec_list = calculate_all_vectors(Template, Match)
    
    if args.rms:
        rms_val = angle_rms(temp_vec_list, match_vec_list)
        print(rms_val)
    elif args.avg:
        mean_val = angle_mean(temp_vec_list, match_vec_list)
        print(mean_val)
    elif args.adt:
        angle_rmsd = angle_difference_test(temp_vec_list, match_vec_list)
        print(angle_rmsd)
    elif args.wadt:
        template_residue_positions = get_residue_positions(Template)
        w_angle_rmsd = weighted_angle_difference_test(temp_vec_list, match_vec_list, template_residue_positions)
        print(w_angle_rmsd)
    else:
        print('Please choose to calculate either angle rms (-rms), angle difference test (-adt) or weighted angle difference test (-wadt)')
