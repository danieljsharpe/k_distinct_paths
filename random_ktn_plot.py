'''
Python script to plot a random kinetic transition network embedded in a 2D space with
energies according to a model potential with noise, dumped from random_ktn.cpp

Daniel J. Sharpe
Feb 2019
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from copy import deepcopy

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
def plot_data(min_ens,min_pos,ts_conns,plot_type):

    xi, yi = np.linspace(-2.,2.,100), np.linspace(-2.,2.,100)
    xi, yi = np.meshgrid(xi,yi)
#    min_posX,min_posY = np.meshgrid(min_pos[:,0],min_pos[:,1])
    rbfi = interpolate.Rbf(min_pos[:,0],min_pos[:,1],min_ens,function="gaussian")
    approx_ens = rbfi(xi,yi)
    plt.imshow(approx_ens,cmap="coolwarm",vmin=-4.,vmax=4.,extent=[xi.min(),xi.max(), \
               yi.min(),yi.max()], origin="lower")
#    if plot_type==0: # plot all points, not joined
    plt.scatter(min_pos[:,0],min_pos[:,1],c=min_ens,cmap="coolwarm",vmin=-4.,vmax=4.)
    if plot_type==1: # plot points pairwise such that connected points are joined
        for conn in ts_conns:
            pt1, pt2 = min_pos[conn[0]-1,:], min_pos[conn[1]-1,:]
#            print "conn[0]-1:", conn[0]-1, "conn[1]-1:", conn[1]-1
#            print "pt1:", pt1, "pt2:", pt2
#            print "d:", ((pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2)**0.5
            plt.plot([pt1[0],pt2[0]],[pt1[1],pt2[1]],"k-")
#            break
    plt.xlim(-2.,2.), plt.ylim(-2.,2.)
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    min_ens = read_data("min.data.dummy",1)
    ts_ens, ts_conns = read_data("ts.data.dummy",2)
    min_pos = read_pos("min.pos.dummy")
    ts_pos = read_pos("ts.pos.dummy")
#    min_ens = [x for (y,x) in sorted(zip(min_ens,
    plot_data(min_ens,min_pos,ts_conns,1)
