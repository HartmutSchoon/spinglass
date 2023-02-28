import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

DEBUG = True
overlap_column_name = 'linked_overlap'
output_filename = './data/test_overlaps'

"""f = open("T_list.txt", "r")
T_list = []
for line in f:
    T_list.append(float(line))

if DEBUG:
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
overlap.to_csv(output_filename, index=False)"""
overlap = pd.read_csv(output_filename)

#Make Histogram
q_start = 0
q_end = 1
n_bins = 50

input_array = overlap["T0.5"].to_numpy()

delta_bin = (q_end-q_start)/n_bins
bin_array = np.array([q_start+delta_bin*x+delta_bin/2 for x in range(n_bins)])

histo, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo = histo/np.sum(histo)
delta_bin = bin_edges[1]-bin_edges[0]
bins = [e+delta_bin/2 for e in bin_edges[:-1]]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(bins, histo, 'kv')
ax.plot(bins, histo, 'k')
ax.set_xlabel("q")
#ax.legend(loc='lower right')
#ax.set_xscale('log')
ax.grid()
fig.savefig(output_filename+'_histo.png')
plt.show() 


print("BP")
    