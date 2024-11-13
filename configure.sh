#!/bin/bash

# Set a predefined GEOCLIM configuration, provided that the files exist in templates
# Practically, this script some current files (config files and some Fortran
# source files) by predefined ones, stored in local template/ directory
#
# 1 input argument: name of configuration 
# Existing predefined configurations are:
#   ref   (ERA5 reanalysis forcings)
#   GFDL  (GFDL climate outputs forcings)
#   IPSL  (IPSL-CM5A2 climate outputs forcings)


cd config/
cp -f templates/IO_CONDITIONS.$1 IO_CONDITIONS
cp -f templates/cond_p20.${1}.dat cond_p20.dat
cp -f templates/GCM_input_conditions.$1 GCM_input_conditions
cd ../source/
cp -f templates/constante.f90.$1 constante.f90
cp -f templates/shape.inc.$1 shape.inc

