#!/bin/bash
#PBS -A UCBK0027
#PBS -N G-car
#PBS -o geoclim-car.log
#PBS -j oe
#PBS -q regular
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=1:ompthreads=1

### ??
export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

### Run OpenMP program
cd ../executable/
./geoclim-car.exe 0 1 3 0 0

# Recursive resubmission
cd ../job/
./submit_chain.sh CONTINUE_RUN
