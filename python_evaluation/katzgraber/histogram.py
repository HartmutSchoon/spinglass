import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


overlap_L4 = pd.read_csv('/home/hatti/Data/spinglass/spinglass_L4/T_sorted/all_overlaps.csv')
#WARNING: DEBUG PATH!
overlap_L6 = pd.read_csv('/home/hatti/Data/spinglass/spinglass_L6/T_sorted/all_overlaps.csv')
overlap_L8 = pd.read_csv('/home/hatti/Data/spinglass/spinglass_L8/T_sorted/all_overlaps.csv')

#Make Histogram
q_start = -1
q_end = 1
n_bins = 100

#delta_bin = (q_end-q_start)/n_bins
#bin_array = np.array([q_start+delta_bin*x+delta_bin/2 for x in range(n_bins)])

input_array = overlap_L4["T0.5"][0:10_000].to_numpy()
histo_L4_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L4_T05= histo_L4_T05/np.sum(histo_L4_T05)

input_array = overlap_L6["T0.5"][0:10_000].to_numpy()
histo_L6_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L6_T05= histo_L6_T05/np.sum(histo_L6_T05)

input_array = overlap_L8["T0.5"][0:10_000].to_numpy()
histo_L8_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L8_T05= histo_L8_T05/np.sum(histo_L8_T05)

delta_bin = bin_edges[1]-bin_edges[0]
bins = [e+delta_bin/2 for e in bin_edges[:-1]]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(bins, histo_L4_T05 , 'kv', label = r'L=4' )
ax.plot(bins, histo_L4_T05, 'k')
ax.plot(bins, histo_L6_T05 , 'ko', label = r'L=6')
ax.plot(bins, histo_L6_T05, 'k')
ax.plot(bins, histo_L8_T05 , 'ks', label = r'L=8')
ax.plot(bins, histo_L8_T05, 'k')


ax.set_xlabel("$\it{q}$")
ax.set_ylabel("$\it{P(q)}$")
ax.legend(loc='upper right')
#ax.set_xscale('log')
ax.grid()
fig.savefig('histo_oneRun.png')
#plt.show() 

input_array = overlap_L4["T0.5"].to_numpy()
histo_L4_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L4_T05= histo_L4_T05/np.sum(histo_L4_T05)

input_array = overlap_L6["T0.5"].to_numpy()
histo_L6_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L6_T05= histo_L6_T05/np.sum(histo_L6_T05)

input_array = overlap_L8["T0.5"].to_numpy()
histo_L8_T05, bin_edges = np.histogram(input_array, bins=n_bins, range=(q_start,q_end))
histo_L8_T05= histo_L8_T05/np.sum(histo_L8_T05)

delta_bin = bin_edges[1]-bin_edges[0]
bins = [e+delta_bin/2 for e in bin_edges[:-1]]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(bins, histo_L4_T05 , 'kv', label = r'L=4')
ax.plot(bins, histo_L4_T05, 'k')
ax.plot(bins, histo_L6_T05 , 'ko', label = r'L=6')
ax.plot(bins, histo_L6_T05, 'k')
ax.plot(bins, histo_L8_T05 , 'ks', label = r'L=8')
ax.plot(bins, histo_L8_T05, 'k')


ax.set_xlabel("$\it{q}$")
ax.set_ylabel("$\it{P(q)}$")
ax.legend(loc='upper right')
#ax.set_xscale('log')
ax.grid()
fig.savefig('histo_all.png')
#plt.show() 



f = open("T_list.txt", "r")
T_vec = []
for line in f:
    T_vec.append(round(float(line),4))
T_vec=np.array(T_vec)

mean_q_L4 = np.zeros(shape=T_vec.shape)
mean_q_L6 = np.zeros(shape=T_vec.shape)
mean_q_L8 = np.zeros(shape=T_vec.shape)

for idx,T in enumerate(T_vec):
    mean_q_L4[idx]=overlap_L4['T'+str(T)].to_numpy().mean()
    mean_q_L6[idx]=overlap_L6['T'+str(T)].to_numpy().mean()
    mean_q_L8[idx]=overlap_L8['T'+str(T)].to_numpy().mean()


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(T_vec, mean_q_L4 , 'kv', label = r'L=4')
ax.plot(T_vec, mean_q_L4, 'k')
ax.plot(T_vec, mean_q_L6 , 'ko', label = r'L=6')
ax.plot(T_vec, mean_q_L6, 'k')
ax.plot(T_vec, mean_q_L8 , 'ks', label = r'L=8')
ax.plot(T_vec, mean_q_L8, 'k')

ax.set_xlabel("$\it{T}$")
ax.set_ylabel("$\it{|q|}$")
ax.legend(loc='upper right')
ax.grid()
fig.savefig('qMean_over_T.png')
plt.show()
print("BP")
    