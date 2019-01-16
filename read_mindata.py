'''
Script to read in a weighted graph (KTN) from min.data and ts.data files
Also contains functions to write Epath.<PATH>, min.data.fastest.<PATH> and min.data.fastest.all files
Daniel J. Sharpe djs244
Jul 2018
'''

import numpy as np

k_B = 1.9872036*10**(-3) # Boltzmann constant / kcal K^{-1} mol^{-1}

def read_data(fhandle,min_or_ts,read_frqs):
    energies = []
    if read_frqs: frqs = []
    if min_or_ts == "TS": connections = []
    for line in fhandle:
        line = line.split()
        energies.append(float(line[0]))
        if read_frqs: frqs.append(float(line[1]))
        if min_or_ts == "TS":
            connections.append((int(line[3]),int(line[4])))
    if min_or_ts == "MIN" and not frqs:
        return np.array(energies,dtype=float)
    elif min_or_ts == "MIN" and frqs:
        return np.array(energies,dtype=float), np.array(frqs,dtype=float)
    elif min_or_ts == "TS" and not frqs:
        return np.array(energies,dtype=float), np.array(connections,dtype=int)
    elif min_or_ts == "TS" and frqs:
        return np.array(energies,dtype=float), np.array(connections,dtype=int), \
               np.array(frqs,dtype=float)

def get_data(read_frqs=False):
    with open("min.data.regrouped","r") as md_f:
        if not read_frqs: min_energies = read_data(md_f,"MIN",read_frqs)
        else: min_energies, min_frqs = read_data(md_f,"MIN",read_frqs)
    with open("ts.data.regrouped","r") as tsd_f:
        if not read_frqs: ts_energies, ts_conns = read_data(tsd_f,"TS",read_frqs)
        else: ts_energies, ts_conns, ts_frqs = read_data(tsd_f,"TS",read_frqs)
    if not read_frqs: return min_energies, ts_energies, ts_conns
    else: return min_energies, ts_energies, ts_conns, min_frqs, ts_frqs

def write_mindatafastest(path, path_no):
    idcs_on_path = sorted([step[0] for step in path])
    with open("min.data.fastest."+str(path_no),"w") as mdf_f:
        for idx in idcs_on_path: mdf_f.write(str(idx)+"\n")

def write_mindatafastestall(in_mdf):
    with open("min.data.fastest.all","w") as mdfa_f:
        for v in in_mdf.iterkeys():
            if in_mdf[v] is True: mdfa_f.write(str(v)+"\n")

def write_epath(path, G, min_energies, ts_energies, path_no):
    with open("Epath."+str(path_no),"w") as epath_f:
        for n_step in range(len(path)):
            # min
            epath_f.write(str((2*n_step)+1)+"\t")
            epath_f.write(str(min_energies[path[n_step][0]-1])+"\t")
            epath_f.write(str(path[n_step][0])+"\n")
            try:
                next_step_min = path[n_step+1][0]
            except IndexError: # reached end of path
                break
            # ts
            epath_f.write(str((2*n_step)+2)+"\t")
            ts_idx = G[path[n_step][0]][next_step_min][1]
            epath_f.write(str(ts_energies[ts_idx-1])+"\t")
            epath_f.write(str(ts_idx)+"\n")

'''Calculate canonical rate constants (Wales adjacency matrix)'''
def calc_rate_const(E_TS,E_min1,E_min2,frq_min1,frq_min2,frq_ts,Natoms,T,s,t, \
                    min1_idx,min2_idx):
#    if (min1_idx != s and min1_idx != t) and (min2_idx != s and min2_idx != t):
    if True:
        ts_cost1 = np.log((np.longdouble(frq_min1)**((3*Natoms)-6) / \
                           np.longdouble(frq_ts)**((3*Natoms)-7))* \
                   np.exp(-(np.longdouble(E_TS)-np.longdouble(E_min1)) / \
                           (np.longdouble(k_B)*T)))
        ts_cost2 = np.log((np.longdouble(frq_min2)**((3*Natoms-6) / \
                           np.longdouble(frq_ts)**(3*Natoms)-7))* \
                   np.exp(-(np.longdouble(E_TS)-np.longdouble(E_min2)) / \
                           (np.longdouble(k_B)*T)))
    else:
        ts_cost1, ts_cost2 = 1., 2.
    return np.float64(ts_cost1), np.float64(ts_cost2)
