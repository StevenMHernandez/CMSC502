#!/usr/bin/env bash

make clean

make mpi

mpirun -np 36 ./mpi
