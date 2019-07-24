'''
Python script to plot a random kinetic transition network embedded in a 2D space with
energies according to a model potential with noise, dumped from random_ktn.cpp

Daniel J. Sharpe
Feb 2019
'''

import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import interpolate
from copy import deepcopy

# for tex
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

''' read positions stored in file min.pos.dummy and ts.pos.dummy '''
def read_pos(posf):

    pos = []
    with open(posf,"r") as pf:
        for line in pf.readlines():
            pos.append([float(x) for x in line.split()])
    return np.array(pos,dtype=float)

''' read energy and connectivity data stored in min.data.dummy and ts.data.dummy '''
def read_data(dataf,sp_type):

    ens = []
    if sp_type==2: conns = []
    with open(dataf,"r") as df:
        for line in df.readlines():
            ens.append(float(line.split()[0]))
            if sp_type==2:
                conns.append([int(line.split()[3]),int(line.split()[4])])
    if sp_type==1:
        return np.array(ens,dtype=float)
    elif sp_type==2:
        return np.array(ens,dtype=float), np.array(conns,dtype=int)


''' make a 2D plot of the data '''
def plot_data(min_ens,min_pos,ts_conns,plot_type,endpt1,endpt2,npath_strt=0,npath_end=0):

    xi, yi = np.linspace(-2.,2.,100), np.linspace(-2.,2.,100)
    xi, yi = np.meshgrid(xi,yi)
#    min_posX,min_posY = np.meshgrid(min_pos[:,0],min_pos[:,1])
    rbfi = interpolate.Rbf(min_pos[:,0],min_pos[:,1],min_ens,function="gaussian")
    approx_ens = rbfi(xi,yi)
    plt.rcParams["xtick.direction"] = "out"
    plt.rcParams["ytick.direction"] = "out"
    plt.imshow(approx_ens,cmap="coolwarm",vmin=-4.,vmax=4.,extent=[xi.min(),xi.max(), \
               yi.min(),yi.max()], origin="lower")
    # plot vertices
    plt.scatter(min_pos[:,0],min_pos[:,1],c=min_ens,cmap="coolwarm",vmin=-4.,vmax=4.)
#    plt.scatter(min_pos[endpt1-1,0],min_pos[endpt1-1,1],marker="ro")
#    plt.scatter(min_pos[endpt2-1,0],min_pos[endpt2-1,1],marker="ro")
    plt.plot([min_pos[endpt1-1,0],min_pos[endpt2-1,0]],[min_pos[endpt1-1,1],min_pos[endpt2-1,1]], \
             "ro",markersize=10,zorder=10)
    if plot_type==1: # plot points pairwise such that connected points are joined
        for conn in ts_conns:
            pt1, pt2 = min_pos[conn[0]-1,:], min_pos[conn[1]-1,:]
            plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],"k-")
    elif plot_type==2: # plot only the edges in the rate-limiting cut
        rate_lim_edges = []
        with open("rate_lim_cut.dat","r") as rlcdat:
            for line in rlcdat.readlines():
                rate_lim_edges.append(int(line.split()[0]))
        for ts_id in rate_lim_edges:
            conn = ts_conns[ts_id-1]
            pt1, pt2 = min_pos[conn[0]-1,:], min_pos[conn[1]-1,:]
            plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],"k-")
    elif plot_type==3: # plot edges that appear in Epath.x files for npaths paths
        for i in range(npath_strt,npath_end+1):
            with open("Epath."+str(i),"r") as epathf:
                nstep = 0
                minid1, minid2 = 0, 0
                for line in epathf.readlines():
                    line = line.split()
                    nstep += 1
                    if nstep%2==0: # transition state
                        ts_id = int(line[2]) # transition state
                    elif nstep!=1: # plot the min-sad-min triplet
                        minid2 = int(line[2])
                        conn = ts_conns[ts_id-1]
                        pt1, pt2 = min_pos[conn[0]-1,:], min_pos[conn[1]-1,:]
                        plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],"k-")
                        minid1 = minid2
                    else: # beginning of pathway
                        minid1 = int(line[2])
    # plot
    plt.xlim(-2.,2.), plt.ylim(-2.,2.)
    ax = plt.gca()
    ax.set_xlabel("$x$",fontsize=20)
    ax.set_ylabel("$y$",fontsize=20)
    ax.tick_params(labelsize=14)
    fig = plt.gcf()
#    fig.set_size_inches(8,8) # set figure size
#    cbar = plt.colorbar()
#    cbar.ax.set_ylabel("$V(x,y)$",fontsize=20)
#    cbar.ax.tick_params(labelsize=14)
    plt.savefig("fig.png",format="png")
    plt.savefig("fig.ps",format="ps")
    plt.show()

if __name__ == "__main__":

    # if = 1: plot all edges, if = 2: plot edges listed in "rate_lim_cut.dat" only, 
    # if = 3: plot edges listed in a range of "Epath.x" files (also need to provide npath_strt and npath_end)
    edge_plot_type = int(sys.argv[1])
    endpt1, endpt2 = int(sys.argv[2]), int(sys.argv[3])
    assert edge_plot_type in [1,2,3]
    if edge_plot_type == 3:
        npath_strt, npath_end = int(sys.argv[4]), int(sys.argv[5])

    min_ens = read_data("min.data.dummy",1)
    ts_ens, ts_conns = read_data("ts.data.dummy",2)
    min_pos = read_pos("min.pos.dummy")
    ts_pos = read_pos("ts.pos.dummy")
#    min_ens = [x for (y,x) in sorted(zip(min_ens,
    if edge_plot_type == 3: plot_data(min_ens,min_pos,ts_conns,edge_plot_type, \
        endpt1,endpt2,npath_strt,npath_end)
    else: plot_data(min_ens,min_pos,ts_conns,edge_plot_type,endpt1,endpt2)
