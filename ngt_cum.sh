#!/bin/bash

# loop over min/ts.data.fastest.x files and calculate NGT rate constant at each step
# usage:
# ./ngt_cum.sh <min_pathno> <max_pathno> <intvl_pathno>

# make sure "NGT 0 F" appears in the pathdata file (and no other calculation keywords)

# min/ts.data.fastest.x files (containing single-col list of min/ts indices, respectively)
#   need to be present in the working directory
# the full database must be stored in min/ts.data.all files, and the endpoint nodes for the full
#   database must be stored in min.A/B.all files (in the working dir)
# need to set a PATHSAMPLE executable <ps_exec> that should output a single-line file "ngt_rateconsts.dat"
# need to set a (Python) script <parse_exec> to write min.data and ts.data files from
#   min/ts.data.fastest.x files at each iteration
# output is written to file "ngt_cum.dat"

min_pathno=$1
intvl_pathno=$2
max_pathno=$3
ps_exec='/home/djs244/PATHSAMPLE_DUMPINFOMAP2'
parse_exec='/home/djs244/UTILITIES/make_datafastest_files.py'

rm PS_ngt.out ngt_rateconsts.dat ngt_cum.dat min.data.fastest.new ts.data.fastest.new # cleanup
for i in $(seq -f "%04g" $min_pathno $intvl_pathno $max_pathno)
do
    echo "Iteration of cumulative NGT calculation: $i"
    python $parse_exec min.data.fastest.$i ts.data.fastest.$i .all
    mv min.data.fastest.new min.data
    mv ts.data.fastest.new ts.data
    mv min.A.new min.A
    mv min.B.new min.B
    $ps_exec > PS_ngt.out
    line=$(sed -n 1p ngt_rateconsts.dat)
    echo "$i    $line" >> "ngt_cum.dat"
    rm ngt_rateconsts.dat
done
