#!/bin/bash
#SBATCH -o zopfli.%u.o%j
#SBATCH -e zopfli.%u.e%j
#SBATCH -n 16
#SBATCH --exclusive
#SBATCH -p normal # normal queue; use development or serial if queue wait time is too long
#SBATCH -t 01:00:00 # 1 hour

for iters in 3 5 15 20;
do
    ./results.sh $iters
done
