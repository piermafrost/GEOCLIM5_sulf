#!/bin/bash

# Set GEOCLIM configuration, provided that the files exist in templates
# 1 input argument: name of configuration 

cp -f templates/IO_CONDITIONS.$1 IO_CONDITIONS
cp -f templates/cond_p20.${1}.dat cond_p20.dat
cp -f templates/GCM_input_conditions.$1 GCM_input_conditions
cd ../source/
cp -f templates/constante.f90.$1 constante.f90
cp -f templates/shape.inc.$1 shape.inc

