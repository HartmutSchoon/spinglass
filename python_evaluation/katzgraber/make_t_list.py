

def calc_T_list():  
    T_list=[]
    T = 0.5
    while T < 2.0:
        T_list.append(T)
        deltaT = 0.035*T + 0.03
        T += deltaT
    return T_list

    
T_list = calc_T_list()
file = open('T_list.txt','w')
for T in T_list:
    file.write(str(T)+"\n")
file.close()