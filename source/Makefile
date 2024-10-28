# This is a makefile for GEOCLIM. Support both intel and GNU compiler.
#
# This makefile includes one optional argument variable: MODE, defining the compilation model
# By default, MODE = standard, but two other modes can be specified while calling make: "debug" and "optim"
#     make MODE=standard    (equivalent to just "make")
#     make MODE=optim       (for fast use)
#     make MODE=debug       (extra debugging options)


#################################################
# Main customization variables (default values) #
#################################################

# Fortran compiler
FC = gfortran

# Compilation mode (default):
MODE = standard

# Path of netCDF library (must contains subrepertories lib/ and include/)
ncpath = /usr


#############################
# Compiler-dependent flags: #
#############################


ifeq ($(findstring pgfortran, $(FC)), pgfortran)

#=========================#
# PGI compiler: pgfortran #
#=========================#

# Fortran language option
lang_flags = -Mfreeform

# Error and warning options:
ifeq ($(MODE), optim)
warn_flags = -Minform=severe
else ifeq ($(MODE), debug)
warn_flags = -Minform=inform
else # assume standard mode
warn_flags = -Minform=warn
endif

# Debugging options:
ifeq ($(MODE), debug)
debug_flags = -g
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -Mbounds
endif
ifeq ($(MODE), debug)
code_flags += -traceback
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0
else ifeq ($(MODE), optim)
optim_flags = -O4 -fast
else
# if none of previous: assume standard mode:
optim_flags = -O2
endif



else ifeq ($(findstring gfortran, $(FC)), gfortran)

#========================#
# GNU compiler: gfortran #
#========================#

# Fortran language options:
lang_flags = -ffree-form -std=gnu

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -Wall -Wconversion-extra -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code
else
warn_flags = -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wconversion -Wunderflow -Wunreachable-code
ifeq ($(MODE), standard)
warn_flags += -Wpedantic -Wunused-dummy-argument
endif
endif

# Debugging options:
ifneq ($(MODE), optim)
debug_flags = -ffpe-trap=invalid,zero,overflow
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -fbounds-check
endif
ifeq ($(MODE), debug)
code_flags += -fbacktrace -g3
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all
else ifeq ($(MODE), optim)
optim_flags = -O3
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



else ifeq ($(findstring ifort, $(FC)), ifort)

#========================
# Intel compiler: ifort #
#========================

# Fortran langauge options:
lang_flags = -free -132

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -warn all
else
warn_flags = -warn general,usage,declaration,truncated_source,ignore_loc
endif

# Debugging options:
#ifeq ($(MODE), debug)
#debug_flags = -debug full
#endif
# Note: this option would be nice, but it sometimes causes compilation crash with:
# "relocation truncated to fit: R_X86_64_32 against symbol `bio_div1_' defined in COMMON section"

# Options at run time
ifeq ($(MODE), debug)
code_flags = -check all -fp-stack-check -traceback
else ifeq ($(MODE), standard)
code_flags = -check bounds
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all -fstack-security-check
else ifeq ($(MODE), optim)
optim_flags = -Ofast
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



endif



########################
# LIBRARY and INCLUDE: #
########################


# NetCDF library:
inc_flags += -I$(ncpath)/include
lib_flags += -L$(ncpath)/lib
netcdf_flags = -lnetcdf -lnetcdff


##############################################
# Complete list of fortran compilation flags #
##############################################

main_flags = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags)
FFLAGS = $(main_flags) $(inc_flags) $(lib_flags) $(netcdf_flags)


#############################################################################

#############
# Root path #
#############

# Check if path assignment uses correct path
# If not, add 'updatepath' to the executable dependencies

# get root path (remove 'source' from current directory path)
pathcommand = character(len=*), parameter:: geoclim_path = "$(CURDIR:source=)"
# get path assignment command (last line of 'path.inc' file)
currpathcommand = $(shell tail -1 path.inc)

ifneq ($(pathcommand), $(currpathcommand))
pathsource = updatepath path.inc
else
pathsource = path.inc
endif


#################################################
# Fortran source files and compilation commands #
#################################################


# source files:
sources = anoxic.f basin_geometry.f bio_frac.f biological_pump.f carb_dep.f carbo.f chimie.f co2_interpolation.f90 constante.f90 cont_weath.f creades.f dc13_speciation.f degassing.f derivs.f diss_oxygen.f DOfA.f dynsoil.f90 dynsoil_create_init_condition.f90 dynsoil_create_output.f90 dynsoil_create_restart.f90 dynsoil_empirical_laws.f90 dynsoil_initialization.f90 dynsoil_lithium.f90 dynsoil_offline_computation.f90 dynsoil_physical_parameters.f90 dynsoil_read_input.f90 dynsoil_steady_state.f90 dynsoil_write_output.f90 eqcte.f facco2.f oceanic_T.f GCM_io_module.f90 geoclim_create_output.f90 geoclim_mainprog.f geoclim_write_output.f90 geographic_create_output.f90 geographic_write_output.f90 get_cont_pixel.f io_module.f90 load_lithology.f lyso.f netcdf_io_module.f90 newton.f newton2.f ocean_atm_flu.f open_ascii_files.f org_dep.f output_netcdf2ascii.f90 pco2_atm.f phbor.f phfunc.f phfunc2.f phosphorite.f Phydrotherm.f physical_units.f90 printf.f read_conditions.f read_deathnote.f90 read_main_IO_conditions.f read_oceanic_temperature.f90 read_veget.f recycle.f rk4.f rtbis.f rtbis2.f salinity.f sea_omega.f strontium_ratio.f terrestrial_dep.f utils.f90 varset.f


# Executable file:
execut = geoclim.exe


# List of objects
objects_0 = $(sources:.f90=.o)
objects = $(objects_0:.f=.o)


# Rule to build executable
${execut}: ${objects}
	$(FC) $^ -o $@ $(FFLAGS)


# Redefine implicit rule to for making objects:
%.o: %.f
	$(FC) $< -o $@ -c $(FFLAGS)
%.o: %.f90
	$(FC) $< -o $@ -c $(FFLAGS)



# File dependencies:
DOfA.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
anoxic.o: constante.o DOfA.o combine_foam.inc coupler.inc shape.inc $(pathsource)
basin_geometry.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
bio_frac.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
utils.o: path.inc
io_module.o: utils.o
biological_pump.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
carb_dep.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
phfunc.o: constante.o eqcte.o combine_foam.inc coupler.inc shape.inc $(pathsource)
newton.o: phfunc.o
rtbis.o: phfunc.o combine_foam.inc coupler.inc shape.inc $(pathsource)
carbo.o: newton.o rtbis.o combine_foam.inc coupler.inc shape.inc $(pathsource)
co2_interpolation.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
dynsoil_physical_parameters.o: shape.inc
dynsoil_empirical_laws.o: dynsoil_physical_parameters.o
dynsoil.o: dynsoil_empirical_laws.o dynsoil_physical_parameters.o
dynsoil_steady_state.o: dynsoil_empirical_laws.o dynsoil_physical_parameters.o
dynsoil_lithium.o: dynsoil_empirical_laws.o dynsoil_physical_parameters.o
facco2.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
cont_weath.o: carbo.o constante.o dynsoil.o dynsoil_lithium.o dynsoil_steady_state.o facco2.o co2_interpolation.o combine_foam.inc coupler.inc shape.inc $(pathsource)
diss_oxygen.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
dc13_speciation.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
degassing.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
oceanic_T.o: co2_interpolation.o constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
ocean_atm_flu.o: chimie.o constante.o eqcte.o phbor.o combine_foam.inc coupler.inc shape.inc $(pathsource)
org_dep.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
pco2_atm.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
phosphorite.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
Phydrotherm.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
recycle.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
salinity.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
sea_omega.o: constante.o lyso.o combine_foam.inc coupler.inc shape.inc $(pathsource)
strontium_ratio.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
creades.o: anoxic.o biological_pump.o bio_frac.o carb_dep.o constante.o cont_weath.o dc13_speciation.o degassing.o diss_oxygen.o oceanic_T.o ocean_atm_flu.o org_dep.o pco2_atm.o phosphorite.o Phydrotherm.o recycle.o salinity.o sea_omega.o strontium_ratio.o combine_foam.inc coupler.inc shape.inc $(pathsource)
varset.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
derivs.o: varset.o combine_foam.inc coupler.inc shape.inc $(pathsource)
get_cont_pixel.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
dynsoil_create_init_condition.o: co2_interpolation.o constante.o dynsoil_empirical_laws.o get_cont_pixel.o combine_foam.inc coupler.inc dynsoil_physical_parameters.o $(pathsource)
dynsoil_create_output.o: io_module.o netcdf_io_module.o coupler.inc
dynsoil_create_restart.o: netcdf_io_module.o coupler.inc dynsoil_physical_parameters.o
dynsoil_initialization.o: io_module.o dynsoil.o
dynsoil_offline_computation.o: dynsoil_steady_state.o coupler.inc dynsoil_physical_parameters.o
dynsoil_read_input.o: dynsoil_create_init_condition.o io_module.o utils.o dynsoil_physical_parameters.o
dynsoil_write_output.o: netcdf_io_module.o
GCM_io_module.o: physical_units.o utils.o combine_foam.inc coupler.inc shape.inc $(pathsource)
geoclim_create_output.o: io_module.o netcdf_io_module.o coupler.inc
geoclim_write_output.o: constante.o netcdf_io_module.o combine_foam.inc coupler.inc shape.inc $(pathsource)
geographic_create_output.o: io_module.o netcdf_io_module.o coupler.inc
geographic_write_output.o: netcdf_io_module.o
load_lithology.o: io_module.o utils.o combine_foam.inc coupler.inc shape.inc $(pathsource)
phfunc2.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
newton2.o: phfunc2.o
open_ascii_files.o: utils.o combine_foam.inc coupler.inc shape.inc $(pathsource)
output_netcdf2ascii.o: constante.o combine_foam.inc coupler.inc shape.inc $(pathsource)
read_deathnote.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
printf.o: dynsoil_create_restart.o dynsoil_offline_computation.o dynsoil_write_output.o geoclim_write_output.o geographic_write_output.o read_deathnote.o combine_foam.inc coupler.inc shape.inc $(pathsource)
read_conditions.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
read_main_IO_conditions.o: utils.o combine_foam.inc coupler.inc shape.inc $(pathsource)
read_oceanic_temperature.o: shape.inc
read_veget.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
rk4.o: derivs.o
rtbis2.o: phfunc2.o combine_foam.inc coupler.inc shape.inc $(pathsource)
terrestrial_dep.o: combine_foam.inc coupler.inc shape.inc $(pathsource)
geoclim_mainprog.o: basin_geometry.o constante.o creades.o dynsoil_create_output.o dynsoil_initialization.o dynsoil_read_input.o GCM_io_module.o geoclim_create_output.o geographic_create_output.o get_cont_pixel.o load_lithology.o open_ascii_files.o output_netcdf2ascii.o printf.o read_conditions.o read_main_IO_conditions.o read_oceanic_temperature.o read_veget.o rk4.o utils.o varset.o combine_foam.inc coupler.inc shape.inc $(pathsource)



#############################################################################


############################
# Update root path command #
############################

# create path.inc file if it doesn't exist
path.inc:
	touch $@

.PHONY: updatepath # make sure to always execute that command when called
updatepath:
	echo '! File automatically generated by `make`' > path.inc
	echo '$(pathcommand)' >> path.inc


####################
# Special commands #
####################

.PHONY: clean clc echo echo_FC

echo:
	@echo "$(FC) $(FFLAGS)"

echo_FC:
	@echo $(FC)

clean:
	rm -f *.o *.mod *__gen_mod* *__genmod*
clc:
	rm -f *.o *.mod *__gen_mod* *__genmod* $(execut)
