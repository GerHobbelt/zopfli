#!/bin/bash

# Submits many Slurm jobs to measure this zopfli thing


for threads in 1 2 4 8;
do
    for rn in 1 2;
    do
        for rd in 20 100;
        do
            # Set build parameters and build a new zopfli version
            export CFLAGS="-DAM_OMP_THREAD_NUM=$threads -DAM_OMP_T_RAND_NUM=$rn -DAM_OMP_T_RAND_DENOM=$rd"
            export CXXFLAGS="-DAM_OMP_THREAD_NUM=$threads -DAM_OMP_T_RAND_NUM=$rn -DAM_OMP_T_RAND_DENOM=$rd"
            make clean
            make
            LLsub result_up.sh
            # just a guess at how long it will take
            sleep 60
        done
    done
done
