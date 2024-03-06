import os

out_file_list = '8_res_templates.txt'                                                                                                                                       
if os.path.exists(out_file_list):
    os.remove(out_file_list)
    
cwd = os.getcwd()

for root, dirs, files in os.walk(cwd, topdown=False):
    for name in files:
        if name[-4:] == '.pdb':
            with open(out_file_list, 'a') as f:
                f.write(os.path.join(cwd, root, name) + '\n')
