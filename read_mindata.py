'''
Script to read in a weighted graph (KTN) from min.data and ts.data files
Also contains functions to write Epath.<PATH>, min.data.fastest.<PATH> and min.data.fastest.all files
djs244
Jul 2018
'''

import numpy as np

def read_data(fhandle,min_or_ts):
    energies = []
    if min_or_ts == "TS": connections = []
    for line in fhandle:
        line = line.split()
        energies.append(float(line[0]))
        if min_or_ts == "TS":
            connections.append((int(line[3]),int(line[4])))
    if min_or_ts == "MIN": return np.array(energies,dtype=float)
    elif min_or_ts == "TS": return np.array(energies,dtype=float), np.array(connections,dtype=int)

def get_data():
    with open("min.data","r") as md_f:
        min_energies = read_data(md_f,"MIN")
    with open("ts.data","r") as tsd_f:
        ts_energies, ts_conns = read_data(tsd_f,"TS")
    return min_energies, ts_energies, ts_conns

def write_mindatafastest():
    pass

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
            try:
                ts_G_idx = [i[0] for i in G[path[n_step][0]]].index(next_step_min)
                ts_idx = G[path[n_step][0]][ts_G_idx][2]
            except ValueError:
                ts_G_idx = [i[0] for i in G[next_step_min]].index(path[n_step][0])
                ts_idx = G[next_step_min][ts_G_idx][2]
            epath_f.write(str(ts_energies[ts_idx-1])+"\t")
            epath_f.write(str(ts_idx)+"\n")
