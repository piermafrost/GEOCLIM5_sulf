# This is a makefile for GEOCLIM. Support both intel and GNU compiler.
# 
# This makefile includes one optional argument variable: MODE, defining the compilation model
# By default, MODE = standard, but two other modes can be specified while calling make: "debug" and "optim"
#     make MODE=standard    (equivalent to just "make")
#     make MODE=optim       (for fast use)
#     make MODE=debug       (extra debugging options)



# Compililation mode:
MODE ?= standard # Default mode


# Macros:
##################################
FC ?= gfortran# default compiler #
##################################
DUMMYFC = $(FC)
# If you directly do `ifeq ($(FC), gfortran)', for some osbcure reasons, it doesn't work. Whereas `ifeq ($(DUMMYFC), gfortran)' works.


# compiler-dependent flags:

                            #========================#
ifeq ($(DUMMYFC), gfortran) # GNU compiler: gfortran #
                            #========================#

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

                         #========================
ifeq ($(DUMMYFC), ifort) # Intel compiler: ifort #
                         #========================

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



########################
# LIBRARY and INCLUDE: #
########################


# NetCDF library:
ncpath = /usr
inc_flags += -I$(ncpath)/include
lib_flags += -L$(ncpath)/lib
netcdf_flags = -lnetcdf -lnetcdff



FFLAGS = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags) $(inc_flags) $(lib_flags) $(netcdf_flags)



#############################################################################


create_uniform_litho: create_uniform_litho.f90
	$(FC) $(FFLAGS) $^ -o $@

create_slope: create_slope.f90
	$(FC) $(FFLAGS) $^ -o $@

create_flat_slope: create_flat_slope.f90
	$(FC) $(FFLAGS) $^ -o $@

reduce_slope: reduce_slope.f90
	$(FC) $(FFLAGS) $^ -o $@

remove_null_slope: remove_null_slope.f90
	$(FC) $(FFLAGS) $^ -o $@

yvesdat2netcdf: yvesdat2netcdf.f90
	$(FC) $(FFLAGS) $^ -o $@

#############################################################################


#dynsoil_empirical_laws.o: dynsoil_empirical_laws.f90
#	$(FC) -c $(FFLAGS) $< -o $@
#
#create_init_cond_equilibrated.o: create_init_cond_equilibrated.f90
#	$(FC) -c $(FFLAGS) $< -o $@


#############################################################################

# File dependencies:

#create_init_cond_equilibrated.o: dynsoil_empirical_laws.o


#############################################################################


check:
	@echo "$(FC) $(FFLAGS)"

.PHONY: clean clc

clean:
	rm -f *.o *.mod *__gen_mod*
clc:
	rm -f *.o *.mod *__gen_mod* create_uniform_litho create_slope create_flat_slope reduce_slope remove_null_slope yvesdat2netcdf
