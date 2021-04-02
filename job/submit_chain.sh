#!/bin/bash

# This is a template bash script for submitting a GEOCLIM run 
# in several times, one after the other, and modifying the solver
# and printing time steps before launching the next run.
# The code need to be compiled and configured, the present script
# will edit the configuration files automatically.
# The successive runs will be renamed: *_1, *_2, *_3, ...
#
# NOTE: if the argument "CONTINUE_RUN" is passed to the present script,
# it will expect the internal configuration to be already done,
# and will edit the GEOCLIM configuration file and launch the next run.
# This trick is used for the recursive submission of GEOCLIM run.
#
# The actual "job" file that launch geoclime run is a separate file,
# stated by the variable JOB_FILE
# This file must be configured for the user's cluster.


#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#
#<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>#

# <><><><><><><><><> #
#  User's parameters #
# <><><><><><><><><> #

# NOTE: the present file, the job file and the run configuration are stored in directory specific
#       to the run launched. It can be freely modified after submitting the job without any
#       consequence on the run


# GEOCLIM's main configuration file (RELATIVE TO GEOCLIM ROOT DIRECTORY!)
GEOCLIM_IO_FILE='config/IO_CONDITIONS'

# GEOCLIM's run name (if undefined, keep the one in main GEOCLIM IO file)
RUN_NAME='.jobtest'

# Job file (configured for the cluster) and submission commmand
JOB_FILE='run_geoclim.sh'
SUBMIT_COMMAND=qsub

# "real" executable file (if undefined, default is ../executable/geoclim.exe)
# (RELATIVE TO GEOCLIM ROOT DIRECTORY!)
EXECUTABLE='executable/geoclim.exe'

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
STOP_TIMES="1 10 100"   # times to stop and relaunch GEOCLIM
COMBINE_DT=""                  # solver time-step for combine 
CONTWTH_NSKIP="1 1 1" # asynchronous coupling of continental weathering module (x time_step)
DYNSOIL_NSKIP="1 1 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
COMBINE_PRINT_NSKIP="1 1 1"        # printing COMBINE outputs every ... x time_step
GEOGRAP_PRINT_NSKIP="2 10 100"     # printing geographic outputs every ... x time_step
DYNSOIL_PRINT_NSKIP="4 20 100"  # printing Dynsoil outputs every ... x time_step
#
## "realistic example"
# STOP_TIMES="1d+5 1d+6 50d+6"   # times to stop and relaunch GEOCLIM
# COMBINE_DT=""                  # solver time-step for combine 
# CONTWTH_NSKIP="100 1000 10000" # asynchronous coupling of continental weathering module (x time_step)
# DYNSOIL_NSKIP="4 2 1"          # asynchronous coupling of DynSoil module (x contwth_nskip x time_step)
# COMBINE_PRINT_NSKIP="4000 40000 400000"        # printing COMBINE outputs every ... x time_step
# GEOGRAP_PRINT_NSKIP="40000 400000 4000000"     # printing geographic outputs every ... x time_step
# DYNSOIL_PRINT_NSKIP="200000 2000000 40000000"  # printing Dynsoil outputs every ... x time_step

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




if [[ $1 != "CONTINUE_RUN" ]] # => Internal configuration
then                          #%%%%%%%%%%%%%%%%%%%%%%%%%%


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
    tail -n $((GEOCLIM_IO_FILE_NLINES - RUN_NAME_LINENUM)) ../$GEOCLIM_IO_FILE >> toto
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

    cp $0        $RUN_NAME/run/
    cp $JOB_FILE $RUN_NAME/run/
    chmod u+x $RUN_NAME/run/*



    # Submit GEOCLIM job
    # ------------------

    # Go in the local "run" directory
    # + + + + + + + #
    cd $RUN_NAME/run/
    # + + + + + + + #

    # Link to "real" executable (job file must run the link "geoclim.exe")
    test "${EXECUTABLE:0:1}" == "/" || EXECUTABLE=../../../$EXECUTABLE
    ln -s -f $EXECUTABLE geoclim.exe

    #=======================#
    $SUBMIT_COMMAND $JOB_FILE
    #./run_geoclim.sh
    #=======================#

    if [ $? -ne 0 ] # if submission failed, put back original GEOCLIM config files
    then
        cd ../../
	test -e .backup/IO_CONDITIONS || mv -f .backup/IO_CONDITIONS ../$GEOCLIM_IO_FILE
	test -e .backup/cond_p20.dat  || mv -f .backup/cond_p20.dat  $CONFIG_FILE
    fi




else # configuration already done => continue run
     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

	#=======================#
	$SUBMIT_COMMAND $JOB_FILE
	#./run_geoclim.sh
	#=======================#

	if [ $? -ne 0 ] # if submission failed, put back original GEOCLIM config files
	then
	    test -e ../../.backup/IO_CONDITIONS || mv -f ../../.backup/IO_CONDITIONS $GEOCLIM_IO_FILE
	    test -e ../../.backup/cond_p20.dat  || mv -f ../../.backup/cond_p20.dat  $CONFIG_FILE
	fi



    else # End of runs, put back original GEOCLIM config files

	test -e ../../.backup/IO_CONDITIONS || mv -f ../../.backup/IO_CONDITIONS $GEOCLIM_IO_FILE
	test -e ../../.backup/cond_p20.dat  || mv -f ../../.backup/cond_p20.dat  $CONFIG_FILE

    fi




fi
