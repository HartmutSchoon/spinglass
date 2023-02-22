import numpy as np
import pandas as pd
from pathlib import Path
import re
import matplotlib.pyplot as plt

def main():
    #all_results_path = '/users/student/xese4803/Data/spinglass/testData'
    all_results_path = '/home/hatti/Data/spinglass'
    all_results_path = Path(all_results_path)
    filtered_data = []
    for results_path in all_results_path.iterdir():
        calc_one_dir(results_path)
        filtered_data.append(pd.read_csv(results_path / 'T05_filtered_data.csv'))
    #print("BP")
      
    """filtered_data = pd.read_csv('./T05_filtered_data.csv')    """

    """ plt.plot(av_data["run"]/1000, av_data["av_linked Overlapp"] )
    plt.plot(av_data["run"]/1000, av_data["av_katz_energy"] )
    plt.xlabel("Sweep")
    plt.grid()
    plt.show() """
    
    plt.plot(filtered_data[0]["run"]/1000, filtered_data[0]["av_linked_overlapp"] )
    plt.plot(filtered_data[0]["run"]/1000, filtered_data[0]["av_katz_energy"] )
    plt.xlabel("Sweep")
    plt.xscale('log')
    plt.grid()
    plt.show()

    #print("hello World")
    

def calc_one_dir(results_path):
     
        file_regex = re.compile('grid.*tsv')
        file_paths =[file for file in results_path.iterdir() if file_regex.match(file.name)]
        
        print("Read results directory "+results_path.name)
        data = [pd.read_csv(file,sep='\t') for file in file_paths]
        
        data = katz_energy_for_list_of_df(data, z=6, J_square=1, num_particles=1000)
        
        filtered_data = filter_by_T_and_average(data,target_T=0.5) 
        filtered_data = append_df(filtered_data, mov_average_over_half(filtered_data['linked_overlapp']))
        filtered_data = append_df(filtered_data, mov_average_over_half(filtered_data['katz_energy']))
    
        filtered_data.to_csv(results_path / 'T05_filtered_data.csv')
    
    

def filter_by_T_and_average(input_data, target_T):
    '''this function produces a resulting DataFrame with the same size
    and column names as the input DataFrames, where each row is the average
    of the input_data rows where the rows temperature is equal to target_T'''
    
    print("Compute T filtered averaged Dataframe")
    resulting_df = pd.DataFrame(columns=input_data[0].columns)
    #Loop over every row
    for row_idx in range(len(input_data[0])):
        print("Compute T filtered averaged Dataframe:"+str(row_idx+1)+"/"+str(len(input_data[0])))
        av_row = pd.DataFrame(0, index=[0], columns= input_data[0].columns)
        num_matching_rows = 0
        for df in input_data:
            row = df.loc[row_idx,:]
            if row['T'] == target_T:
                av_row += row
                num_matching_rows += 1
        av_row /= num_matching_rows
        resulting_df = pd.concat([resulting_df, av_row], axis=0)
    return resulting_df

            
    
def average_list_of_df(input_data):
    '''Takes a list of pd.Dataframes with equal columns and size and computes
        a DataFrame with columns and size equal to the input_data containing the
        average all DataFrames in input_data.
        resulting_DF[column_name][idx] = mean(df[column_name][idx]) over all df in input_data '''

    result_df = pd.DataFrame(0, index=input_data[0].index, columns=input_data[0].columns)
    for df in input_data:
        result_df += df
    return result_df/len(input_data)
        
 
def mov_average_for_list_of_df(data, column_name):
    '''Loops over a list of of DataFrames and computes the 'mov_average_over_half()' for the specified 
    column in 'column_name' for every DataFrame'''
    
    for (idx,df) in enumerate(data):
        series = df[column_name]
        av_series = mov_average_over_half(series)
        new_df = append_df(df, av_series)
        data[idx] = new_df
    return data

def katz_energy_for_list_of_df(data, z, J_square, num_particles):
    '''Loops over a list of of DataFrames and computes the 'katz_energy_for_df' for the specified 
    column in 'column_name' for every DataFrame'''
    print("Compute katz energy for all data")
    for (idx,df) in enumerate(data):
        print("Compute katz energy for all data:" + str(idx+1) + "/" + str(len(data)))
        katz_series = katz_energy_for_df(df, z, J_square, num_particles)
        new_df = append_df(df,katz_series)
        data[idx] = new_df
    return data


def append_df(input_df, input_series):
    '''function to concat a series to DataFrame'''
    input_df.reset_index(drop=True, inplace=True)
    input_series.reset_index(drop=True, inplace=True)
    return pd.concat([input_df, input_series], axis=1)

def katz_energy_for_df(input_df, z, J_square, num_particles):
    '''adds a column to input df containing the 'katzgraber energy' with column name 'katz_energy' '''
    katz_energy_list = []
    for (idx,row) in input_df.iterrows():
        T = row['T']
        U = row['energy']
        katz_energy_list.append(katz_energy(T,U,z,J_square, num_particles))
    katz_energy_series = pd.Series(katz_energy_list, name= 'katz_energy')
    return katz_energy_series
    
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
        
def katz_energy(T, U, z, J_square, num_particles):
    '''computes the 'Energy' Term from the Katzgraber Paper
    our Energy ist normalized yet-> we need to do it here'''
    return 1-2.0*T*abs(U)/(z*J_square*num_particles)
        
if __name__ == "__main__":
    main()