import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

DEBUG = False
overlap_column_name = 'overlap'
path_to_data = '/home/hatti/Data/spinglass/spinglass_L8/T_sorted'
output_file = Path(path_to_data) / 'all_overlaps.csv'

f = open("T_list.txt", "r")
T_list = []
for line in f:
    T_list.append(round(float(line),4))

if DEBUG:
    T_list=T_list[0:3]
    
data_all_T = []    
for T in T_list:
    #Load Data
    results_path = Path(path_to_data)
    results_path = results_path / ('T' + str(T))
    file_paths = [file for file in results_path.iterdir()]
    data = [pd.read_csv(path).drop('Unnamed: 0', axis=1) for path in file_paths]


    #there are some dataframes which weren't stored correctly. Therefore they get dropped!
    bad_idx = []
    for idx, df in enumerate(data):
        if len(df)!= 10000:
            bad_idx.append(idx)
    #List needs to be selected in reverse, so the indexing in data wont be disturbed
    for idx in sorted(bad_idx, reverse=True):
        del data[idx]
        
    #write down all
    data_all_T.append((T,data))

#get all overlap values for each T in one big DF
#overlap_df = pd.DataFrame()

n_overlap_values_per_df = len(data_all_T[0][1][0])

pd_columns = []
overlap_mat = np.zeros((0))
for (T,data) in data_all_T:
    print("T: " + str(T))
    overlap_this_T = np.zeros((0,1))
    for df in data:
        overlap_this_T=np.append(overlap_this_T, df[overlap_column_name])
    series=pd.Series(data = overlap_this_T, name="T"+str(T))
    pd_columns.append(series)


overlap = pd.concat(pd_columns, axis=1, keys=[s.name for s in pd_columns])
overlap.to_csv(output_file, index=False)
print("BP")