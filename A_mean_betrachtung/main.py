import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    data = pd.read_table("T_linear_0_4.0.tsv")
    plt.plot(data["T"], data["A_mean"] )
    plt.xlabel("T")
    plt.ylabel("<A>")
    plt.show()
    
    """ sim_delta_t=t_test()
    plt.plot(sim_delta_t["T"], sim_delta_t["delta_T"] )
    plt.xlabel("T")
    plt.ylabel("dT")
    plt.show()
    print(sim_delta_t) """
    
    """ logistic_define() """
    
    print("BP")
    
def logistic_define():
    T_start = 0.5
    T_end = 4
    T = np.linspace(T_start, T_end, num= 100)
    dT = 100/(1+np.exp(4.0*(T-1.75)))
    df_dict = {
        "T" : T,
        "dT": dT
    }
    data = pd.DataFrame(df_dict)
    
    plt.plot(data["T"], data["dT"])
    plt.xlabel("T")
    plt.ylabel("dT")
    plt.show()
    print("BP")
    
def t_test():
    K = 100
    t_start = 0.5
    t_end = 2
    beta_start = 1/t_start
    beta_end = 1/t_end
    
    delta_beta = (beta_end-beta_start)/(K-1.0)
    
    current_beta = beta_start
    
    t_list = []
    delta_t_list = []
    beta_list = []
    delta_beta_list = []
    
    for i in range(K):
        current_T = 1/current_beta         
        t_list.append(current_T)
        beta_list.append(current_beta)
        current_beta += delta_beta
        
    #t_list.reverse()
    for idx in  range(1,len(t_list)):
        delta_t_list.append(t_list[idx]-t_list[idx-1])
        delta_beta_list.append(beta_list[idx]-beta_list[idx-1])
    delta_t_list.append(delta_t_list[-1])  
    delta_beta_list.append(delta_beta_list[-1]) 
    
    df_dict = {
        'T':t_list,
        'delta_T':delta_t_list,
        'Beta': beta_list,
        'delta_Beta':delta_beta_list}
    sim_delta_t= pd.DataFrame(df_dict)
    return sim_delta_t
        
    

if __name__ == "__main__":
    main()
