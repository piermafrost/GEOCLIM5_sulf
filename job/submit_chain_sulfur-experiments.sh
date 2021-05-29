#!/bin/bash

# This is a template bash script for submitting a GEOCLIM run 
# in several times, one after the other, and modifying the solver
# and printing time steps before launching the next run.
# The code need to be compiled and configured, the present script
# will edit the configuration files automatically.
# The successive runs will be renamed: *_1, *_2, *_3, ...
#
# The actual "job" file that launch geoclime run is a separate file,
# stated by the variable JOB_FILE (usually, run_geoclim.sh)
# This file must be configured for the user's cluster, if submitted as
# a batch process. Templates can be found in the current repertory.
#
# This script use the log file .config-queue (in main GEOCLIM directory)
# to resolve conflict of access to GEOCLIM configuration files.
# You can type `./submit_chain.sh ERASE_CONFIG_QUEUE` to deleted it
# in case an error occur and the status is "busy" while no run is waiting
# to start.
#
# If you want to restore the original configuration files (normally done at
# the end of ALL sucessive runs launched by this script), type
# `./submit_chain.sh RESTORE_CONFIG`.


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

# <><><><><><><><><> #
#  User's parameters #
# <><><><><><><><><> #

# NOTE: the present file, the job file and the run configuration are stored in directory specific
#       to the run launched. It can be freely modified after submitting the job without any
#       consequence on the run



# Sulfur experiment run set-ups (on Cheyenne cluster)
# <><><><><><><><><><><><><><><><><><><><><><><><><><>


# Run-independent variables
# -------------------------

# GEOCLIM's main configuration file (RELATIVE TO GEOCLIM ROOT DIRECTORY!)
GEOCLIM_IO_FILE='config/IO_CONDITIONS'

# Job file (configured for the cluster) and submission commmand
# If you want to simply run the model, and not submitting a batch process, set JOB_FILE='run_geoclim_basic.sh' and SUBMIT_COMMAND=''
JOB_FILE='run_geoclim.sh'
SUBMIT_COMMAND='qsub' # note: '' (empty variable) means direct execution of the job file, and not submitting it as a batch process 

LOG_FILE="geoclim.log"

# Dynsoil initial conditions
# 4 possibilities:
#    * 'startup:eq' (start from automatically-computed equilibrium condition)
#    * 'startup:null' (start from null regolith)
#    * name_of_init_file (start from given restart file)
#    * '' (empty variable: keep restart mode and file from GEOCLIM main config file)
DYNSOIL_INIT='startup:eq'


# Run-dependent variable:
# -----------------------

# Abrupt perturbation (3 successive runs)
STOP_TIMES="1d+5 1d+6 50d+6"   # times to stop and relaunch GEOCLIM
COMBINE_DT="0.25 0.25 0.25"    # solver time-step for combine 
CONTWTH_NSKIP="100 2000 40000" # asynchronous coupling of continental weathering module (x time_step)
DYNSOIL_NSKIP="4 2 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
COMBINE_PRINT_NSKIP="4000   40000   400000"    # printing COMBINE outputs every ... x time_step
GEOGRAP_PRINT_NSKIP="40000  400000  4000000"   # printing geographic outputs every ... x time_step
DYNSOIL_PRINT_NSKIP="200000 2000000 40000000"  # printing Dynsoil outputs every ... x time_step

# FOR CARB TRADE-OFF (cannot have asynchronous DynSoil coupling because carbonate wth is computed outside DynSoil)
#DYNSOIL_NSKIP="1 1 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)

# Progressive perturbation (1 single run)
#STOP_TIMES="50d+6"    # times to stop and relaunch GEOCLIM
#COMBINE_DT="0.25"     # solver time-step for combine 
#CONTWTH_NSKIP="40000" # asynchronous coupling of continental weathering module (x time_step)
#DYNSOIL_NSKIP="1"     # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
#COMBINE_PRINT_NSKIP="400000"    # printing COMBINE outputs every ... x time_step
#GEOGRAP_PRINT_NSKIP="4000000"   # printing geographic outputs every ... x time_step
#DYNSOIL_PRINT_NSKIP="40000000"  # printing Dynsoil outputs every ... x time_step


###  Carb sulf wth - O2 fdbk: P-Hyd++, land++
#COMBINE_INIT='restart/geoclim/output.calib_Phyd-land-fdbk-plus'
#
#RUN_NAME='.PyrW+50-Carb_Phyd-land-fdbk-plus'
#EXECUTABLE='executable/geoclim-car_Phydlandfdbkplus.exe'
#RUN_NAME='.PyrW+50-Carb-prog_Phyd-land-fdbk-plus'
#EXECUTABLE='executable/geoclim-car-prog_Phydlandfdbkplus.exe'

#    Sil sulf wth - O2 fdbk: P-Hyd++, land++
#RUN_NAME='.PyrW+50-Sil_Phyd-land-fdbk-plus'
#EXECUTABLE='executable/geoclim-sil_Phydlandfdbkplus.exe'
#RUN_NAME='.PyrW+50-Sil-prog_Phyd-land-fdbk-plus'
#EXECUTABLE='executable/geoclim-sil-prog_Phydlandfdbkplus.exe'

###  Carb sulf wth - O2 fdbk: P-Hyd++, land+
#COMBINE_INIT='restart/geoclim/output.calib_Phyd-land-fdbk'
#
#RUN_NAME='.PyrW+50-Carb_Phyd-land-fdbk'
#EXECUTABLE='executable/geoclim-car_Phydlandfdbk.exe'
#RUN_NAME='.PyrW+50-Carb-prog_Phyd-land-fdbk'
#EXECUTABLE='executable/geoclim-car-prog_Phydlandfdbk.exe'

#    Sil sulf wth - O2 fdbk: P-Hyd++, land+
#RUN_NAME='.PyrW+50-Sil_Phyd-land-fdbk'
#EXECUTABLE='executable/geoclim-sil_Phydlandfdbk.exe'
#RUN_NAME='.PyrW+50-Sil-prog_Phyd-land-fdbk'
#EXECUTABLE='executable/geoclim-sil-prog_Phydlandfdbk.exe'

###  Carb sulf wth - O2 fdbk: P-Hyd++
#COMBINE_INIT='restart/geoclim/output.calib_Phyd-fdbk-plus'
#
#RUN_NAME='.PyrW+50-Carb_Phyd-fdbk-plus'
#EXECUTABLE='executable/geoclim-car_Phydfdbkplus.exe'
#RUN_NAME='.PyrW+50-Carb-prog_Phyd-fdbk-plus'
#EXECUTABLE='executable/geoclim-car-prog_Phydfdbkplus.exe'

#    Sil sulf wth - O2 fdbk: P-Hyd++
#RUN_NAME='.PyrW+50-Sil_Phyd-fdbk-plus'
#EXECUTABLE='executable/geoclim-sil_Phydfdbkplus.exe'
#RUN_NAME='.PyrW+50-Sil-prog_Phyd-fdbk-plus'
#EXECUTABLE='executable/geoclim-sil-prog_Phydfdbkplus.exe'

###  Carb sulf wth - O2 fdbk: P-Hyd+
#COMBINE_INIT='restart/geoclim/output.calib_Phyd-fdbk'
#
#RUN_NAME='.PyrW+50-Carb_Phyd-fdbk'
#EXECUTABLE='executable/geoclim-car_Phydfdbk.exe'
#RUN_NAME='.PyrW+50-Carb-prog_Phyd-fdbk'
#EXECUTABLE='executable/geoclim-car-prog_Phydfdbk.exe'

#    Sil sulf wth - O2 fdbk: P-Hyd+
#RUN_NAME='.PyrW+50-Sil_Phyd-fdbk'
#EXECUTABLE='executable/geoclim-sil_Phydfdbk.exe'
#RUN_NAME='.PyrW+50-Sil-prog_Phyd-fdbk'
#EXECUTABLE='executable/geoclim-sil-prog_Phydfdbk.exe'

###  Carb sulf wth - O2 fdbk: no P
#COMBINE_INIT='restart/geoclim/output.calib-noPfdbk'
#
#RUN_NAME='.PyrW+50-Carb_noPfdbk'
#EXECUTABLE='executable/geoclim-car_noPfdbk.exe'
#RUN_NAME='.PyrW+50-Carb-prog_noPfdbk'
#EXECUTABLE='executable/geoclim-car-prog_noPfdbk.exe'

#    Sil sulf wth - O2 fdbk: no P
#RUN_NAME='.PyrW+50-Sil_noPfdbk'
#EXECUTABLE='executable/geoclim-sil_noPfdbk.exe'
#RUN_NAME='.PyrW+50-Sil-prog_noPfdbk'
#EXECUTABLE='executable/geoclim-sil-prog_noPfdbk.exe'

###  Carb sulf wth - O2 fdbk: no P, sed ML-
#COMBINE_INIT='restart/geoclim/output.calib-noPfdbk-red'
#
#RUN_NAME='.PyrW+50-Carb_noPfdbk-red'
#EXECUTABLE='executable/geoclim-car_noPfdbkred.exe'
#RUN_NAME='.PyrW+50-Carb-prog_noPfdbk-red'
#EXECUTABLE='executable/geoclim-car-prog_noPfdbkred.exe'

#    Sil sulf wth - O2 fdbk: no P, sed ML-
#RUN_NAME='.PyrW+50-Sil_noPfdbk-red'
#EXECUTABLE='executable/geoclim-sil_noPfdbkred.exe'
#RUN_NAME='.PyrW+50-Sil-prog_noPfdbk-red'
#EXECUTABLE='executable/geoclim-sil-prog_noPfdbkred.exe'

###  Carb sulf wth - no O2 fdbk
#COMBINE_INIT='restart/geoclim/output.GFDL_PI_eq'
#
#RUN_NAME='.PyrW+50-Carb_noO2fdbk'
#EXECUTABLE='executable/geoclim-car_noO2fdbk.exe'
#RUN_NAME='.PyrW+50-Carb-prog_noO2fdbk'
#EXECUTABLE='executable/geoclim-car-prog_noO2fdbk.exe'

#    Sil sulf wth - no O2 fdbk
#RUN_NAME='.PyrW+50-Sil_noO2fdbk'
#EXECUTABLE='executable/geoclim-sil_noO2fdbk.exe'
#RUN_NAME='.PyrW+50-Sil-prog_noO2fdbk'
#EXECUTABLE='executable/geoclim-sil-prog_noO2fdbk.exe'

###  RESTART with standard O2 feedback
COMBINE_INIT='restart/geoclim/output.GFDL_PI_eq'

###  Carb sulf wth - standard O2 feedback
RUN_NAME='.PyrW+50-Carb'
EXECUTABLE='executable/geoclim-car.exe'
#RUN_NAME='.PyrW+50-Carb-prog'
#EXECUTABLE='executable/geoclim-car-prog.exe'

#    Sil sulf wth - standard O2 feedback
#COMBINE_INIT='restart/geoclim/output.GFDL_PI_eq'
#RUN_NAME='.PyrW+50-Sil'
#EXECUTABLE='executable/geoclim-sil.exe'
#RUN_NAME='.PyrW+50-Sil-prog'
#EXECUTABLE='executable/geoclim-sil-prog.exe'

###  Carb sulf wth - Pyrite & Kerogen perturbation, without add. P (standard O2 feedback)
#RUN_NAME='.PyrKerW+10-noP-Carb'
#EXECUTABLE='executable/geoclim-car-ker+10-noP.exe'
#RUN_NAME='.PyrKerW+10-noP-Carb-prog'
#EXECUTABLE='executable/geoclim-car-ker+10-noP-prog.exe'

###  Carb sulf wth - Pyrite & Kerogen perturbation, with add. P (standard O2 feedback)
#RUN_NAME='.PyrKerW+10-P-Carb'
#EXECUTABLE='executable/geoclim-car-ker+10-P.exe'
#RUN_NAME='.PyrKerW+10-P-Carb-prog'
#EXECUTABLE='executable/geoclim-car-ker+10-P-prog.exe'

###  Carb sulf wth - trade-off
#RUN_NAME='.PyrW+50-Carb-trdf'
#EXECUTABLE='executable/geoclim-car-trdf.exe'

#    Sil sulf wth - trade-off
#RUN_NAME='.PyrW+50-Sil-trdf'
#EXECUTABLE='executable/geoclim-sil-trdf.exe'

#    H2SO4 leaching
#RUN_NAME='.PyrW+50-H2SO4'
#EXECUTABLE='executable/geoclim-h2s.exe'

###  Carb sulf wth - fixed oceanic temperature
#RUN_NAME='.PyrW+50-Carb_cst-oceT'
#EXECUTABLE='executable/geoclim-car-oceTcste.exe'

###  Carb sulf wth - fixed CO2
#!! NOTE: YOU NEED TO UPDATE config/GCM_io_condition TO SPECIFY 1 CO2 LEVEL
#RUN_NAME='.PyrW+50-Carb_cst-CO2'
#EXECUTABLE='executable/geoclim-car-CO2cste.exe'


# <><><><><><><><><><><><><><><><><><><><><><><><><><>

# Solver and printing parameters
# must be *same length* "list" of values (lenght = number of successive runs)
# Undefined of empty variable => keep the values already in the configuration file

#STOP_TIMES="1d+5 1d+6 50d+6"   # times to stop and relaunch GEOCLIM
#COMBINE_DT="0.25 0.25 0.25"    # solver time-step for combine 
#CONTWTH_NSKIP="100 2000 40000" # asynchronous coupling of continental weathering module (x time_step)
#DYNSOIL_NSKIP="4 2 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
#COMBINE_PRINT_NSKIP="4000   40000   400000"    # printing COMBINE outputs every ... x time_step
#GEOGRAP_PRINT_NSKIP="40000  400000  4000000"   # printing geographic outputs every ... x time_step
#DYNSOIL_PRINT_NSKIP="200000 2000000 40000000"  # printing Dynsoil outputs every ... x time_step

### EXAMPLES:
#
## basic test
#STOP_TIMES="1 10 100"   # times to stop and relaunch GEOCLIM
#COMBINE_DT=""                  # solver time-step for combine 
#CONTWTH_NSKIP="1 1 1" # asynchronous coupling of continental weathering module (x time_step)
#DYNSOIL_NSKIP="1 1 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
#COMBINE_PRINT_NSKIP="1 1 1"        # printing COMBINE outputs every ... x time_step
#GEOGRAP_PRINT_NSKIP="2 10 100"     # printing geographic outputs every ... x time_step
#DYNSOIL_PRINT_NSKIP="4 20 100"  # printing Dynsoil outputs every ... x time_step

#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#



## FUNCTIONS
############

function number_of_lines() {
    cat $1 | wc -l
}

function get_second_arg() {
    out=$2
    out=${out//\"/}
    out=${out//\'/}
    echo $out
}

function is_commented() {
    line=$1
    line=${line// }                    # remove all blanks
    line=${line%%\#*}                  # remove longest comment
    test -z $line && echo 0 || echo 1  # if remaining string is empty => commented line => return 0
}

function read_lines() {
    # syntax: read_lines file line_start [line_stop]
    # by default, line_stop = line_start
    file=$1
    start=$2
    test -z $3 && stop=$start || stop=$3
    count=$((stop - start + 1))
    head -n $stop "$file" | tail -n $count
}





## MAIN PROGRAM
###############



# Select case according to argument passed to the script, or presence of "CONTINUE_RUN"
# =====================================================================================



if [ "$1" == "ERASE_CONFIG_QUEUE" ]
then
    rm -f ../.config-queue
    test $? -eq 0 && echo "Warning: queue file for configuration access deleted."
    exit 0
fi



#======================================================================#



if [ ! -e CONTINUE_RUN ]
then
    # No file "CONTINUE_RUN" in the current repertory => Internal configuration
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    # Read GEOCLIM main IO file
    # get run name, 2nd config file name, output directory, restarts directory and restarts
    # -------------------------------------------------------------------------------------

    ls ../$GEOCLIM_IO_FILE > /dev/null || exit 1


    nline=0
    nuncomment=0
    GEOCLIM_IO_FILE_NLINES=`number_of_lines ../$GEOCLIM_IO_FILE`

    # line #1 > run name 
    while [ $nuncomment -lt 1 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    test -z $RUN_NAME && RUN_NAME=`get_second_arg $line`
    RUN_NAME_LINENUM=$nline

    # line #2 > output directory
    while [ $nuncomment -lt 2 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    OUTPUTDIR=`get_second_arg $line`

    # line #3 > 2nd config file ("cond_p20.dat")
    while [ $nuncomment -lt 3 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    CONFIG_FILE=`get_second_arg $line`

    # line #5 > current COMBINE restart
    while [ $nuncomment -lt 5 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    file=`get_second_arg $line`
    COMBINE_RESTARTDIR="${file%\/*}/"
    COMBINE_INIT_LINENUM=$nline

    # line #34 > COMBINE restart "root" name
    while [ $nuncomment -lt 34 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    COMBINE_RESTART=`get_second_arg $line`

    # line #47 > DynSoil restart "root" name
    while [ $nuncomment -lt 47 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    DYNSOIL_RESTART=`get_second_arg $line`

    # line #48 > DynSoil restart mode
    while [ $nuncomment -lt 48 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    DYNSOIL_RESTART_MODE_LINENUM=$nline

    # line #49 > current DynSoil restart
    while [ $nuncomment -lt 49 ]
    do
        nline=$((nline + 1))
        line=`read_lines ../$GEOCLIM_IO_FILE $nline`
	nuncomment=$((nuncomment + `is_commented "$line"`))
    done
    line=${line//\'}
    line=${line//\"}
    DYNSOIL_RESTARTDIR="${line%\/*}/"
    DYNSOIL_INIT_LINENUM=$nline



    # Signal to restore configuration file
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if [ "$1" == "RESTORE_CONFIG" ]
    then
	mv -f .backup/IO_CONDITIONS ../$GEOCLIM_IO_FILE
	mv -f .backup/cond_p20.dat  ../$CONFIG_FILE
	exit 0
    fi



    # Check that a GEOCLIM run is not currently trying to access the configuration files
    # ----------------------------------------------------------------------------------

    if [ -f ../.config-queue ]
    then
        status_string=`tail -n 1 ../.config-queue`
	while [ "$status_string" == "busy" ]
	do
	    echo "configuration files busy. Run set-up postponed for 1 minute"
	    sleep 1m
            status_string=`tail -n 1 ../.config-queue`
	done
    fi

    # Signal that the configuration files are now busy for a run set-up
    echo "busy" >> ../.config-queue



    # Create local storage directory
    # ------------------------------

    rm -r -f $RUN_NAME
    mkdir $RUN_NAME
    mkdir $RUN_NAME/config
    mkdir $RUN_NAME/run
    test -d .backup || mkdir .backup # backup directory common to all runs



    # Write GEOCLIM IO info in local configuration storage directory
    # --------------------------------------------------------------

    ls ../$CONFIG_FILE > /dev/null || exit 1

    echo $GEOCLIM_IO_FILE > $RUN_NAME/config/geoclim_io_file
    echo $CONFIG_FILE > $RUN_NAME/config/config_file

    echo $GEOCLIM_IO_FILE_NLINES       >> $RUN_NAME/config/IOfile_structure
    echo $RUN_NAME_LINENUM             >> $RUN_NAME/config/IOfile_structure
    echo $COMBINE_INIT_LINENUM         >> $RUN_NAME/config/IOfile_structure
    echo $DYNSOIL_RESTART_MODE_LINENUM >> $RUN_NAME/config/IOfile_structure
    echo $DYNSOIL_INIT_LINENUM         >> $RUN_NAME/config/IOfile_structure

    echo $OUTPUTDIR >> $RUN_NAME/config/outputdir

    echo $COMBINE_RESTARTDIR >> $RUN_NAME/config/combine_restart
    echo $COMBINE_RESTART    >> $RUN_NAME/config/combine_restart

    echo $DYNSOIL_RESTARTDIR >> $RUN_NAME/config/dynsoil_restart
    echo $DYNSOIL_RESTART    >> $RUN_NAME/config/dynsoil_restart

    echo $EXECUTABLE > $RUN_NAME/config/executable

    echo $RUN_NAME > $RUN_NAME/config/run_name

    echo 1 > $RUN_NAME/config/run_number


    for k in $STOP_TIMES
    do
        echo $k >> $RUN_NAME/config/stop_times
    done
    test -e $RUN_NAME/config/stop_times && tend=`head -n 1 $RUN_NAME/config/stop_times`

    for k in $COMBINE_DT
    do
        echo $k >> $RUN_NAME/config/combine_dt
    done
    test -e $RUN_NAME/config/combine_dt && dt=`head -n 1 $RUN_NAME/config/combine_dt`

    for k in $CONTWTH_NSKIP
    do
        echo $k >> $RUN_NAME/config/contwth_nskip
    done
    test -e $RUN_NAME/config/contwth_nskip && dt_contwth=`head -n 1 $RUN_NAME/config/contwth_nskip`

    for k in $DYNSOIL_NSKIP
    do
        echo $k >> $RUN_NAME/config/dynsoil_nskip
    done
    test -e $RUN_NAME/config/dynsoil_nskip && dt_dynsoil=`head -n 1 $RUN_NAME/config/dynsoil_nskip`

    for k in $COMBINE_PRINT_NSKIP
    do
        echo $k >> $RUN_NAME/config/combine_print_nskip
    done
    test -e $RUN_NAME/config/combine_print_nskip && tprint=`head -n 1 $RUN_NAME/config/combine_print_nskip`

    for k in $GEOGRAP_PRINT_NSKIP
    do
        echo $k >> $RUN_NAME/config/geograp_print_nskip
    done
    test -e $RUN_NAME/config/geograp_print_nskip && tprint_geog=`head -n 1 $RUN_NAME/config/geograp_print_nskip`

    for k in $DYNSOIL_PRINT_NSKIP
    do
        echo $k >> $RUN_NAME/config/dynsoil_print_nskip
    done
    test -e $RUN_NAME/config/dynsoil_print_nskip && tprint_dynsoil=`head -n 1 $RUN_NAME/config/dynsoil_print_nskip`


    # Determine path relative to current directory
    test "${CONFIG_FILE:0:1}" == "/" || CONFIG_FILE=../$CONFIG_FILE



    # Backup of GEOCLIM configuration files
    # -------------------------------------

    test -e .backup/IO_CONDITIONS || cp ../$GEOCLIM_IO_FILE .backup/IO_CONDITIONS
    test -e .backup/cond_p20.dat  || cp $CONFIG_FILE        .backup/cond_p20.dat



    # Re-create GEOCLIM main IO file, with modified run name (add "_1")
    # -----------------------------------------------------------------

    rm -f toto
    read_lines ../$GEOCLIM_IO_FILE 1 $((RUN_NAME_LINENUM - 1)) > toto
    echo "RUN_NAME: '${RUN_NAME}_1'" >> toto
    if [ -z "$COMBINE_INIT" ]
    then
        read_lines ../$GEOCLIM_IO_FILE $((RUN_NAME_LINENUM + 1)) $((DYNSOIL_RESTART_MODE_LINENUM - 1))     >> toto
    else
        read_lines ../$GEOCLIM_IO_FILE $((RUN_NAME_LINENUM + 1)) $((COMBINE_INIT_LINENUM - 1))             >> toto
        echo "COMBINE_INIT_FILE: '$COMBINE_INIT'"                                                          >> toto
        read_lines ../$GEOCLIM_IO_FILE $((COMBINE_INIT_LINENUM + 1)) $((DYNSOIL_RESTART_MODE_LINENUM - 1)) >> toto
    fi
    case "$DYNSOIL_INIT" in
        "") # empty variable => keep cofig file as it is
            tail -n $((GEOCLIM_IO_FILE_NLINES - DYNSOIL_RESTART_MODE_LINENUM + 1)) ../$GEOCLIM_IO_FILE         >> toto
            ;;
        "startup:eq"|"startup:null") # startup modes, do not need to edit restart file name
            echo "'$DYNSOIL_INIT'"                                                                             >> toto
            tail -n $((GEOCLIM_IO_FILE_NLINES - DYNSOIL_RESTART_MODE_LINENUM)) ../$GEOCLIM_IO_FILE             >> toto
            ;;
        *) # None of above => consider variable as init file name
            echo "'restart'"                                                                                   >> toto
            read_lines ../$GEOCLIM_IO_FILE $((DYNSOIL_RESTART_MODE_LINENUM + 1)) $((DYNSOIL_INIT_LINENUM - 1)) >> toto
            echo "'$DYNSOIL_INIT'"                                                                             >> toto
            tail -n $((GEOCLIM_IO_FILE_NLINES - DYNSOIL_INIT_LINENUM)) ../$GEOCLIM_IO_FILE                     >> toto
            ;;
    esac
    #============================#
    mv -f toto ../$GEOCLIM_IO_FILE
    #============================#



    # Re-create 2nd config file (cond_p20) with current timesteps
    # -----------------------------------------------------------

    rm -f toto
    read_lines $CONFIG_FILE 1 19  >  toto
    test -z $dt             && read_lines $CONFIG_FILE 20 >> toto || echo $dt             >> toto
    echo 0d0                      >> toto
    test -z $tend           && read_lines $CONFIG_FILE 22 >> toto || echo $tend           >> toto
    test -z $tprint         && read_lines $CONFIG_FILE 23 >> toto || echo $tprint         >> toto
    test -z $tprint_geog    && read_lines $CONFIG_FILE 24 >> toto || echo $tprint_geog    >> toto
    echo -1.                      >> toto
    read_lines $CONFIG_FILE 26 28 >> toto
    test -z $dt_contwth     && read_lines $CONFIG_FILE 29 >> toto || echo $dt_contwth     >> toto
    test -z $dt_dynsoil     && read_lines $CONFIG_FILE 30 >> toto || echo $dt_dynsoil     >> toto
    test -z $tprint_dynsoil && read_lines $CONFIG_FILE 31 >> toto || echo $tprint_dynsoil >> toto
    tail -n 3 $CONFIG_FILE        >> toto
    #=====================#
    mv -f toto $CONFIG_FILE
    #=====================#



    # Make a copy of present file and job file in the local storage directory
    # -----------------------------------------------------------------------

    cp $0        $RUN_NAME/run/submit_chain.sh
    cp $JOB_FILE $RUN_NAME/run/
    chmod u+x $RUN_NAME/run/*
    # !!Signal to the script's copy that next run will be a resubmission (with configuration already done)!!
    touch $RUN_NAME/run/CONTINUE_RUN



    # Submit GEOCLIM job
    # ------------------

    # Go in the local "run" directory
    # + + + + + + + #
    cd $RUN_NAME/run/
    # + + + + + + + #

    # Link to "real" executable (job file must run the link "geoclim.exe")
    test "${EXECUTABLE:0:1}" == "/" || EXECUTABLE=../../../$EXECUTABLE
    ln -s -f $EXECUTABLE geoclim.exe

    # Modify the log file of the job file (if asked)
    test -z $LOG_FILE || perl -pi -e "s/geoclim.log/$LOG_FILE/g" $JOB_FILE

    #===========================#
    if [ -z $SUBMIT_COMMAND ]
    then
        ./$JOB_FILE
    else
        $SUBMIT_COMMAND $JOB_FILE
    fi
    #===========================#

    if [ $? -ne 0 ] # if submission failed, put back original GEOCLIM config files
    then
        cd ../../
	test -e .backup/IO_CONDITIONS && mv -f .backup/IO_CONDITIONS ../$GEOCLIM_IO_FILE
	test -e .backup/cond_p20.dat  && mv -f .backup/cond_p20.dat  $CONFIG_FILE
    fi



#======================================================================#



else
    # configuration already done => continue run
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #
    # WARNING: This part is supposed to be executed in the copy the present file
    # that is in ./$RUN_NAME/run/


    # Get configuration info
    # (expect temporary storage files to be already written)
    # ------------------------------------------------------

    GEOCLIM_IO_FILE=`cat ../config/geoclim_io_file`
    CONFIG_FILE=`cat ../config/config_file`
    GEOCLIM_IO_FILE_NLINES=`read_lines       ../config/IOfile_structure 1`
    RUN_NAME_LINENUM=`read_lines             ../config/IOfile_structure 2`
    COMBINE_INIT_LINENUM=`read_lines         ../config/IOfile_structure 3`
    DYNSOIL_RESTART_MODE_LINENUM=`read_lines ../config/IOfile_structure 4`
    DYNSOIL_INIT_LINENUM=`read_lines         ../config/IOfile_structure 5`
    OUTPUTDIR=`cat ../config/outputdir`
    COMBINE_RESTARTDIR=`read_lines ../config/combine_restart 1`
    COMBINE_RESTART=`read_lines    ../config/combine_restart 2`
    DYNSOIL_RESTARTDIR=`read_lines ../config/dynsoil_restart 1`
    DYNSOIL_RESTART=`read_lines    ../config/dynsoil_restart 2`
    EXECUTABLE=`cat ../config/executable`
    RUN_NAME=`cat ../config/run_name`
    RUN_NUMBER=`cat ../config/run_number`


    # Get time steps info
    test -e ../config/stop_times && tstart=`head -n 1 ../config/stop_times`

    the_end=0
    rm -f toto

    if [ -e ../config/stop_times ]
    then
        tail -n+2 ../config/stop_times > toto
        mv toto ../config/stop_times
        test -s ../config/stop_times && tend=`head -n 1 ../config/stop_times` || the_end=1
    fi

    if [ -e ../config/combine_dt ]
    then
        tail -n+2 ../config/combine_dt > toto
        mv toto ../config/combine_dt
        test -s ../config/combine_dt && dt=`head -n 1 ../config/combine_dt` || the_end=1
    fi

    if [ -e ../config/contwth_nskip ]
    then
        tail -n+2 ../config/contwth_nskip > toto
        mv toto ../config/contwth_nskip
        test -s ../config/contwth_nskip && dt_contwth=`head -n 1 ../config/contwth_nskip` || the_end=1
    fi

    if [ -e ../config/dynsoil_nskip ]
    then
        tail -n+2 ../config/dynsoil_nskip > toto
        mv toto ../config/dynsoil_nskip
        test -s ../config/dynsoil_nskip && dt_dynsoil=`head -n 1 ../config/dynsoil_nskip` || the_end=1
    fi

    if [ -e ../config/combine_print_nskip ]
    then
        tail -n+2 ../config/combine_print_nskip > toto
        mv toto ../config/combine_print_nskip
        test -s ../config/combine_print_nskip && tprint=`head -n 1 ../config/combine_print_nskip` || the_end=1
    fi

    if [ -e ../config/geograp_print_nskip ]
    then
        tail -n+2 ../config/geograp_print_nskip > toto
        mv toto ../config/geograp_print_nskip
        test -s ../config/geograp_print_nskip && tprint_geog=`head -n 1 ../config/geograp_print_nskip` || the_end=1
    fi

    if [ -e ../config/dynsoil_print_nskip ]
    then
        tail -n+2 ../config/dynsoil_print_nskip > toto
        mv toto ../config/dynsoil_print_nskip
        test -s ../config/dynsoil_print_nskip && tprint_dynsoil=`head -n 1 ../config/dynsoil_print_nskip` || the_end=1
    fi


    # Determine path relative to current directory
    test "${GEOCLIM_IO_FILE:0:1}"    == "/"  ||  GEOCLIM_IO_FILE=../../../$GEOCLIM_IO_FILE
    test "${CONFIG_FILE:0:1}"        == "/"  ||  CONFIG_FILE=../../../$CONFIG_FILE
    test "${OUTPUTDIR:0:1}"          == "/"  ||  OUTPUTDIR=../../../$OUTPUTDIR
    test "${COMBINE_RESTARTDIR:0:1}" == "/"  &&  LOC_COMBINE_RESTARTDIR=$COMBINE_RESTARTDIR  ||  LOC_COMBINE_RESTARTDIR=../../../$COMBINE_RESTARTDIR
    test "${DYNSOIL_RESTARTDIR:0:1}" == "/"  &&  LOC_DYNSOIL_RESTARTDIR=$DYNSOIL_RESTARTDIR  ||  LOC_DYNSOIL_RESTARTDIR=../../../$DYNSOIL_RESTARTDIR



    # Move restart files
    # ------------------

    mv ${OUTPUTDIR}${COMBINE_RESTART}${RUN_NAME}_${RUN_NUMBER}    ${LOC_COMBINE_RESTARTDIR}
    test $? -ne 0 && the_end=1
    mv ${OUTPUTDIR}${DYNSOIL_RESTART}${RUN_NAME}_${RUN_NUMBER}.nc ${LOC_DYNSOIL_RESTARTDIR}
    test $? -ne 0 && the_end=1



    if [ $the_end -eq 0 ]
    then



	# Check that a GEOCLIM run is not currently trying to access the configuration files
	if [ -f ../../../.config-queue ]
	then
	    status_string=`tail -n 1 ../../../.config-queue`
	    while [ "$status_string" == "busy" ]
	    do
		echo "configuration files busy. Run set-up postponed for 1 minute"
		sleep 1m
		status_string=`tail -n 1 ../../../.config-queue`
	    done
	fi

	# Signal that the configuration files are now busy for a run set-up
	echo "busy" >> ../../../.config-queue



	# Re-create GEOCLIM main IO file, with modified run and restart names
	# -------------------------------------------------------------------

	rm -f toto
	read_lines $GEOCLIM_IO_FILE 1 $((RUN_NAME_LINENUM - 1)) > toto
	echo "RUN_NAME: '${RUN_NAME}_$((RUN_NUMBER + 1))'" >> toto
	read_lines $GEOCLIM_IO_FILE $((RUN_NAME_LINENUM + 1)) $((COMBINE_INIT_LINENUM - 1)) >> toto
	echo "COMBINE_INIT_FILE:    '$COMBINE_RESTARTDIR${COMBINE_RESTART}${RUN_NAME}_${RUN_NUMBER}'" >> toto
	read_lines $GEOCLIM_IO_FILE $((COMBINE_INIT_LINENUM + 1)) $((DYNSOIL_RESTART_MODE_LINENUM - 1)) >> toto
	echo "'restart'" >> toto
	read_lines $GEOCLIM_IO_FILE $((DYNSOIL_RESTART_MODE_LINENUM + 1)) $((DYNSOIL_INIT_LINENUM - 1)) >> toto
	echo "'$DYNSOIL_RESTARTDIR${DYNSOIL_RESTART}${RUN_NAME}_${RUN_NUMBER}.nc'" >> toto
	tail -n $((GEOCLIM_IO_FILE_NLINES - DYNSOIL_INIT_LINENUM)) $GEOCLIM_IO_FILE >> toto
	#=========================#
	mv -f toto $GEOCLIM_IO_FILE
	#=========================#



	# Re-create 2nd config file (cond_p20) with current timesteps
	# -----------------------------------------------------------

	rm -f toto
	read_lines $CONFIG_FILE 1 19  >  toto
	test -z $dt             && read_lines $CONFIG_FILE 20 >> toto || echo $dt             >> toto
	test -z $tstart         && read_lines $CONFIG_FILE 22 >> toto || echo $tstart         >> toto
	test -z $tend           && read_lines $CONFIG_FILE 22 >> toto || echo $tend           >> toto
	test -z $tprint         && read_lines $CONFIG_FILE 23 >> toto || echo $tprint         >> toto
	test -z $tprint_geog    && read_lines $CONFIG_FILE 24 >> toto || echo $tprint_geog    >> toto
	echo -1.                      >> toto
	read_lines $CONFIG_FILE 26 28 >> toto
	test -z $dt_contwth     && read_lines $CONFIG_FILE 29 >> toto || echo $dt_contwth     >> toto
	test -z $dt_dynsoil     && read_lines $CONFIG_FILE 30 >> toto || echo $dt_dynsoil     >> toto
	test -z $tprint_dynsoil && read_lines $CONFIG_FILE 31 >> toto || echo $tprint_dynsoil >> toto
	tail -n 3 $CONFIG_FILE        >> toto
	#=====================#
	mv -f toto $CONFIG_FILE
	#=====================#



	# Resubmit GEOCLIM job
	# --------------------

	# Increase run number
	RUN_NUMBER=$((RUN_NUMBER + 1))
	rm -f ../config/run_number
	echo $RUN_NUMBER > ../config/run_number

	#===========================#
	if [ -z $SUBMIT_COMMAND ]
	then
	    ./$JOB_FILE
	else
	    $SUBMIT_COMMAND $JOB_FILE
	fi
	#===========================#

	if [ $? -ne 0 ] # if submission failed, put back original GEOCLIM config files
	then
	    test -e ../../.backup/IO_CONDITIONS && mv -f ../../.backup/IO_CONDITIONS $GEOCLIM_IO_FILE
	    test -e ../../.backup/cond_p20.dat  && mv -f ../../.backup/cond_p20.dat  $CONFIG_FILE
	fi



    else # End of runs, put back original GEOCLIM config files

	test -e ../../.backup/IO_CONDITIONS && mv -f ../../.backup/IO_CONDITIONS $GEOCLIM_IO_FILE
	test -e ../../.backup/cond_p20.dat  && mv -f ../../.backup/cond_p20.dat  $CONFIG_FILE

    fi



fi

