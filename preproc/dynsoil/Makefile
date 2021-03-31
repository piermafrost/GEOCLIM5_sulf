# This is a makefile for fortran using the compiler and netcdf library currently loaded on the Cheyenne cluster.
# Cheyenne environmental variables used: FC, and NETCDF
# 
# This makefile includes one optional argument variable: MODE, defining the compilation model
# By default, MODE = standard, but two other modes can be specified while calling make: "debug" and "optim"
#     make MODE=standard    (equivalent to just "make")
#     make MODE=optim       (for fast use)
#     make MODE=debug       (extra debugging options)
# 
# WARNING: ALL FILES SHOULD HAVE THE SAME EXTENSION (SPECIFIED BY THE VARIABLE extension)




# Compililation mode:
MODE ?= standard


# Macros:
FC ?= gfortran # A priori, already an environmental variable of the cluster


# compiler-dependent flags:

                       ##########################
ifeq ($(FC), gfortran) # GNU compiler: gfortran #
                       ##########################

	# Fortran language options:
	#lang_flags = -ffree-form -std=f95
	lang_flags = -ffree-form -std=f2003 # some syntax of my programs are not legal in fortran 95 standard

	# Error and warning options:
	ifeq ($(MODE), debug)
		warn_flags = -pedantic-errors -Wall -Wconversion-extra -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code
	else
		warn_flags = -Wpedantic -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wtabs -Winteger-division -Wconversion -Wimplicit-interface -Wimplicit-procedure -Wundefined-do-loop -Wunderflow -Wunused-dummy-argument -Wunreachable-code
	endif

	# Debugging options:
	ifneq ($(MODE), optim)
		debug_flags = -ffpe-trap=invalid,zero,overflow,underflow -g3
	endif

	# Code generation options:
	ifneq ($(MODE), optim)
		code_flags = -fbounds-check
	endif

	# Optimization options:
	ifeq ($(MODE), debug)
		optim_flags = -O0 -fstack-protector-all
	endif
	ifeq ($(MODE), optim)
		optim_flags = -O3
	endif
	# if none of previous: assume standard mode:
	optim_flags ?= -O1

endif

                    #########################
ifeq ($(FC), ifort) # Intel compiler: ifort #
                    #########################

	# Fortran langauge options:
	lang_flags = -free -132

	# Error and warning options (these are the ones I think necessary, some of them are already set by default):
	ifeq ($(MODE), debug)
		warn_flags = -warn all
	else
		warn_flags = -warn general,usage,declaration,truncated_source,interfaces,ignore_loc
	endif

	# Debugging options:
	ifeq ($(MODE), debug)
		debug_flags = -debug full
	endif

	# Options at run time
	ifeq ($(MODE), debug)
		code_flags = -check all -fp-stack-check -traceback
	endif

	# Optimization options:
	ifeq ($(MODE), debug)
		optim_flags = -O0 -fstack-protector-all -fstack-security-check
	endif
	ifeq ($(MODE), optim)
		optim_flags = -Ofast
	endif
	# if none of previous: assume standard mode:
	optim_flags ?= -O1

endif



# NetCDF library:
ncpath = /usr#  $(NETCDF) # -> environmental variable of the cluster storing the netcdf path
inc_flags += -I$(ncpath)/include
lib_flags += -L$(ncpath)/lib
netcdf_flags = -lnetcdf -lnetcdff



FFLAGS = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags) $(inc_flags) $(lib_flags) $(netcdf_flags)


#############################################################################


create_init_cond_null: create_init_cond_null.o
	$(FC) $(FFLAGS) $^ -o $@

create_init_cond_equilibrated: create_init_cond_equilibrated.o dynsoil_empirical_laws.o
	$(FC) $(FFLAGS) $^ -o $@

#############################################################################


# Generic object-making rule:
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@


#############################################################################

# File dependencies:

create_init_cond_equilibrated.o: dynsoil_empirical_laws.o updatelink
create_init_cond_null.o: updatelink
dynsoil_empirical_laws.o: updatelink


#############################################################################


.PHONY: clean check updatelink

updatelink:
	./link.sh

check:
	echo "$(FC) $(FFLAGS)"

clean:
	rm -f *.mod *.o *__genmod*
