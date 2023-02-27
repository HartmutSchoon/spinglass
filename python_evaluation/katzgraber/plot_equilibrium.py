from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def mov_average_over_half(input):
    '''takes a pd.Series and calculates for every index in the series
    the average of all elements between the index and half of the index
    The average is then stored in the output Series
    result[idx] = mean(input[idx/2:idx])'''
    
    results = [];
    for idx in range(len(input)):
        idx_half = round(idx/2)
        results.append(input[idx_half:idx+1].mean())
    return pd.Series(results, name='av_'+input.name)

#------------ VARS START ------------
T = 0.5
num_particles = 10*10*10
z = 6
J_squared = 1
#------------ VARS END------------


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
#data = data[40:52]

#Calc moving average over data
for idx,df in enumerate(data):
    print("Moving average for " + str(idx+1)+ "/" + str(len(data)))
    df['av_energy']=mov_average_over_half(df['energy'])
    df['av_linked_overlap']=mov_average_over_half(df['linked_overlap'])

#Calc katzgraber energy
for df in data:
    print("Katzgraber Energy for " + str(idx+1)+ "/" + str(len(data)))
    df['katz_energy'] =1.0 - 2.0*T*abs(df['av_energy'])/(z*J_squared*num_particles);


#calc mean values over all simulation runs
linked_overlap = np.zeros((len(data[0]),len(data)))
katz_energy = np.zeros((len(data[0]),len(data)))
sweep = data[0]["run"].to_numpy()/1000;

for idx, df in enumerate(data):
    linked_overlap[:,idx]=df["av_linked_overlap"]
    katz_energy[:,idx]=df["katz_energy"]
        
mean_katz = np.mean(katz_energy, axis=1)
mean_overlap = np.mean(linked_overlap, axis=1)

std_err_katz = np.std(katz_energy, axis=1) / np.sqrt(np.shape(katz_energy)[1])
std_err_overlap = np.std(linked_overlap, axis=1) / np.sqrt(np.shape(linked_overlap)[1])

np.savez('equilib_data.npz',sweep=sweep,
            mean_katz=mean_katz,mean_overlap=mean_overlap,
            std_err_katz=std_err_katz, std_err_overlap=std_err_overlap)


output =np.load('equilib_data.npz')
sweep = output['sweep']
mean_katz=output['mean_katz']
mean_overlap=output['mean_overlap']
std_err_katz=output['std_err_katz']
std_err_overlap=output['std_err_overlap']

""" fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(data[0]['run']/1000,data[0]['katz_energy'],
            linestyle='solid', color='black', label = r'$1-2T|U|/z$' )
ax.plot(data[0]['run']/1000,data[0]['av_linked_overlap'], 
            linestyle='dashdot', color='black', label = r'$\langle q_l \rangle$')
ax.set_xlabel("Sweep")
ax.legend(loc='lower right')
ax.set_xscale('log')
ax.grid()
fig.savefig('T05_10E6Sweeps.png')
plt.show() """


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(sweep,mean_katz,
            linestyle='solid', color='black', label = r'$1-2T|U|/z$' )
ax.fill_between(sweep, mean_katz+std_err_katz,mean_katz-std_err_katz, color='gray')
ax.plot(sweep,mean_overlap, 
            linestyle='dashdot', color='black', label = r'$\langle q_l \rangle$')
ax.fill_between(sweep, mean_overlap+std_err_overlap,mean_overlap-std_err_overlap, color='gray')
ax.set_xlabel("Sweep")
ax.legend(loc='lower right')
ax.set_xscale('log')
ax.grid()
fig.savefig('T'+str(T)+'_10E6Sweeps.png')
plt.show() 




