#!/bin/bash
#PBS -A UCBK0027
#PBS -N GEOCLIM
#PBS -o geoclim.log
#PBS -j oe
#PBS -q regular
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1:ompthreads=1

### ??
export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run OpenMP program
./geoclim.exe 0 1 3 0 0

# Recursive resubmission
test $? -eq 0 && ./submit_chain.sh

