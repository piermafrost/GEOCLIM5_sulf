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


#################################################
# Fortran source files and compilation commands #
#################################################


# source files:
sources = compare_template.f90 netcdf_io_module.f90 utils.f90

# Executable file:
execut = compare_template.exe


objects = $(sources:.f90=.o)



# Rule to build executable
${execut}: ${objects}
	$(FC) $^ -o $@ $(FFLAGS)


# Redefine implicit rule to for making objects:
%.o: %.f90
	$(FC) $< -o $@ -c $(FFLAGS)



# File dependencies:
compare_template.o: netcdf_io_module.o utils.o path.inc shape.inc



#############################################################################


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
