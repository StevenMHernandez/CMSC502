#!/usr/bin/env bash

#
# Run this file as an executable as `./timing_test.sh &` so that it runs in the background!
#

make clean
make mpi

#echo "\c" > tmp/timing_runs.txt
#for ((p=1;p<=6;p++)); do
#    for ((i=10;i<=25;i++)); do
#        sqr=`expr $p \\* $p`
#        mpirun -np $sqr ./mpi $i >> tmp/timing_runs.txt
#    done
#done

echo "\c" > tmp/accuracy_runs.txt
for ((p=1;p<=6;p++)); do
    for ((i=1;i<=16;i++)); do
        sqr=`expr $p \\* $p`
        mpirun -np $sqr ./mpi $i >> tmp/accuracy_runs.txt
    done
done