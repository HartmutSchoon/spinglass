import numpy as np
import pandas as pd
from pathlib import Path

f = open("T_list.txt", "r")
T_list = []
for line in f:
    T_list.append(float(line))

T_list=T_list[0:2]
data_all_T = []    
for T in T_list:
    #Load Data
    results_path = Path('/home/hatti/Data/spinglass/T_sorted')
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
    data_all_T.append((T,data))


print("BP")
    