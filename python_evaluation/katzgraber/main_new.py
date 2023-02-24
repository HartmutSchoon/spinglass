import numpy as np
import pandas as pd
from pathlib import Path
import sys
import re
import os


num_equalT_grids = 4

#all_results_path = '/users/student/xese4803/Data/spinglass/testData'
all_results_path = '/home/hatti/Data/spinglass'
all_results_path = Path(all_results_path)
regex = re.compile('results.*')
simulation_paths =[path for path in all_results_path.iterdir() if regex.match(path.name)]

T_sorted_path = all_results_path / 'T_sorted';

os.mkdir(T_sorted_path)
T = float(sys.argv[1])
#T_list = calc_T_list()

# for T in T_list:
this_T_path = T_sorted_path/ ('T' + str(T))
os.mkdir(this_T_path)
for sim_idx,path in enumerate(simulation_paths):
    regex = re.compile('grid.*tsv')
    sim_files = [file for file in path.iterdir() if regex.match(file.name)]
    T_id = 1
    data_one_sim = [pd.read_csv(file,sep='\t') for file in sim_files]
    
    num_lines = data_one_sim[0].shape[0]
    
    matching_T_data = []
    for i in range(num_equalT_grids):
        matching_T_data.append(np.zeros((num_lines,3)))
    
    for line_idx in range(num_lines):
        print("Processing T: " + str(T)
            +" Sim Nr: " + str(sim_idx)
            +" Line Nr: " + str(line_idx))
        found_at_T = 0;
        for data in data_one_sim:
            row = data.loc[line_idx,:]
            if row['T'] == T:
                matching_T_data[found_at_T][line_idx,0] = row['run']
                matching_T_data[found_at_T][line_idx,1] = row['energy']
                matching_T_data[found_at_T][line_idx,2] = row['linked_overlapp']
                found_at_T+=1
        if found_at_T != 4:
            raise ValueError('something went wrong with indexing the np tables')
        found_at_T = 0
    
    for idx,data in enumerate(matching_T_data):
        df = pd.DataFrame(data=data, columns = ['run','energy','linked_overlap'])
        save_name = "data"+"_Sim"+str(sim_idx)+"_Grid"+str(idx)
        df.to_csv(this_T_path / save_name)
        
    print("BP")
        
        

def calc_T_list():  
    T_list=[]
    T = 0.5
    while T < 2.0:
        T_list.append(T)
        deltaT = 0.035*T + 0.03
        T += deltaT
    return T_list
