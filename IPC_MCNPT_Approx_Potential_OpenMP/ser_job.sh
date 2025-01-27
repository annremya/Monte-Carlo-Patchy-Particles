#!/bin/sh
export OMP_NUM_THREADS=4
g++ -O3 -fopenmp IPC_NPT_approx_parallel.cpp
time ./a.out 1 1 > output1.out
exit 0
