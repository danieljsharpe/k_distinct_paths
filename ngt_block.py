'''
Python script to drive PATHSAMPLE NGT rate constant calculations, blocking rate-limiting transition states
(determined by the KDISTINCTPATHS algorithm and listed in the file "rate_lim_cut.dat") from the network
one at a time.
Make sure the pathdata input file includes the keyword "NGT 0 F". Use a modified PATHSAMPLE executable that
dumps the NGT rate constants to the file "ngt_rateconsts.dat"

July 2019
Daniel J. Sharpe
'''

import subprocess
import os.path

### set the string for the location of the PATHSAMPLE executable, and other inputs ###
ps_exec = "/home/djs244/PATHSAMPLE_DUMPINFOMAP"
out_fname = "PS_ngt.out" # name of output file for PATHSAMPLE to write to
max_rle = -1 # set to > 0 to impose a limit on the no. of RLEs considered


if not max_rle > 0:
    res = subprocess.check_output(["wc","-l","rate_lim_cut.dat"])
    n_rle = int(res.split()[0])
    print "Number of rate-limiting edges:", n_rle
else:
    n_rle = max_rle

rle_idcs = [0]*n_rle # list of indices for rate-limiting transition states, in order
i = 0
with open("rate_lim_cut.dat","r") as rlc_f:
    for line in rlc_f.readlines():
        rle_idcs[i] = int(line.split()[0])
        i += 1

if os.path.isfile("ngt_block.dat"): raise IOError # check intended output file does not already exist
ngt_outf = open("ngt_block.dat","w")
for i in range(n_rle+1):
    ps_outf = open(out_fname,"w")
    print "Iteration of NGT with blocking:", i
    if i != 0:
        # delete the rate-limiting transition state of the current iteration from the database
        subprocess.call(["sed","-i",str(rle_idcs[i-1])+"d","ts.data"])
        # subsequent rate-limiting TS indices need to be consistent with new line numbering
        for j, rle_idx in enumerate(rle_idcs[i:]):
            if rle_idcs[i+j] > rle_idcs[i-1]: rle_idcs[i+j] -= 1
    # run NGT calculation
    process = subprocess.Popen([ps_exec],stdout=ps_outf)
    process.wait()
    # process output
    ps_outf.close()
    ngt_ks = subprocess.check_output(["sed","-n","1p","ngt_rateconsts.dat"]) # get first line from file
    ngt_outf.write("%4i   %s" % (i, ngt_ks))
    subprocess.call(["rm","ngt_rateconsts.dat"])
    subprocess.call(["rm",out_fname])
ngt_outf.close()
