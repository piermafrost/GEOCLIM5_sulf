# GEOCLIM5.2

## Updates from GEOCLIM5:
* Possibility to read area, temperature and runoff inputs directly from GCM annual climatology files,
  with some automatic unit conversions.
* Possibility to automatically create DynSoil initial conditions (null regolith, or at equilibrium with initial pCO2).
* Can specify uniform lithology, without input netCDF file
* Add checks for axis matching, missing points, invalid runoff and slope, invalid lithology, and physical units.
  Interactively ask the users if error detected (by default).
* Reorganize the main IO file ('config/IO_CONDITION'), can now add commented (#) and blank lines
* Minor code improvements (e.g., Runge-Kutta 4 scheme, biological 13C fractionation formula...)
* Add simplified sulfur cycle (oceanic SO4^2-, sulfide weathering and sulfate-reduction).


## Required software
* Fortran compiler (note: the Makefile is configured for gfortran, ifort or pgfortran)
* netCDF-Fortran library (see section "How to install netCDF-Fortran library")


## Recommended but not required software
* 'make' software (installed by default on all UNIX/Mac OS)
* Python: some scripts for pre-processing (generate boundary conditions) and visualization are in Python.
  Those scripts use the following packages:
    * numpy
    * netCDF4
    * pylab (some plotting script)
    * matplotlib (plotting scripts)
    * cartopy (a few plotting scripts)
* PyFerret


## Templates
A couple of template GEOCLIM runs are defined. They consist of a set of input data and corresponding configuration files.
The script `make_test` does the complete configuration (pre and post compilation), compiles and runs the desired template
(the pre-existing configuration and source files will be restored when the run is complete).
The expected GEOCLIM outputs for those templates are available in "OUTPUT/templates/", so you can check that you obtain
the same results.
Type `./make_test` for more information.


## `build_GEOCLIM` command
`build_GEOCLIM` is a bash script, meant to edit the needed source files for pre-compilation configuration (data shape,
activated modules, ...) and compile the code. The post-compilation configuration is left to do (which file to use as
input data, names of output and solver parameters, see section "After-compilation configuration").
It uses the Makefile to compile the code.

This command supports a large number of options (specify compiler and options to use...).
Type `./build_GEOCLIM --help` to get detailed information on how to use it.


## Pre-compilation configuration
If you do not use the `build_GEOCLIM command`, here is the list of configuration steps you should follow:
* Which modules to use (file 'source/coupler.inc')
* Number of CO2 levels: 'nclimber' parameter in 'source/shape.inc'
* Geographic resolution (must be consistent with input data): 'nlon' and 'nlat' parameters in 'source/shape.inc'
* Number of lithological classes: variable 'nlitho' in 'source/shape.inc'. Traditionally, 6 classes are used.
  The subroutine 'source/cont_weath.f' explicitly assumes that carbonate rocks is the last lithology class (regardless of
  the number of classes).
  If you want to define more (or less) lithology classes, you should:
    * Leave carbonates as last lithology
    * Update the lithology-dependent parameters in 'source/combine_foam.inc' ('rsw_litho', 'CaMg_rock' and 'OC_in_rock')
    * Update the lithology-dependent parameters in 'source/dynsoil_physical_parameters.inc'
    * Specify which lithology corresponds to "basalt" (variable 'BASALT_LITHO_NUM' in 'source/combine_foam.inc')

  Note that even though 'CaMg_rock' and DynSoil parameters are *not used* if GEOCLIM is not coupled with DynSoil,
  they must be defined consistently with 'nlitho' or the compilation will fail. You may simply fill them with zeros.
* DynSoil vertical resolution: parameter 'nDSlev' in 'source/shape.inc' (must be defined even if not coupled with DynSoil).
  Traditionally 10 vertical levels are used, but it can be changed just by modifying that number.
  If not coupled with DynSoil, its value is not used.
* If coupled to biodiversity module: number of groups (primary producers, predators, ...), in "source/combine_foam.inc"
* The number of oceanic basins, as well as the number of geochemical variables, are not meant to be customized.
* If you don't compile the code with the Makefile, you should also check that the file 'source/path.inc' exists and states
  the correct path of GEOCLIM root directory (with the line `character(len=*), parameter:: geoclim_path = "..."`). If you do
  use the Makefile, it will be updated automatically.

The command `build_GEOCLIM` does all that pre-configuration, except the biodiversity module parameters.

All those parameters must be consistent with the initialization and forcing files (initial ocean chemistry, climate files...)

#### Fixed CO2 run
There is a special case in the model configuration. If there is only 1 CO2 level (nclimber=1), it will be run in "fixed CO2 mode".
This means the amount of CO2 in the atmospheric reservoir will be held constant at the value of the unique CO2 level, whatever
the carbon fluxes. The concentration of the various forms of carbon in the other reservoirs will adjust freely to the atmospheric
concentration (ocean-atmosphere diffusion) and to the carbon fluxes.

This is useful for calibration runs where one wants to hold the atmospheric CO2 constant and adjust the degassing to balance the
silicate weathering flux.


## Compilation

The command `build_GEOCLIM` does the pre-configuration and the compilation. However, it may sometime be useful to compile
the code "manually", for instance if you need to recompile with different optimization options.

#### With Makefile:
First, you need to make sure that the pre-compilation configuration is completed (if `build_GEOCLIM` was not invoked).
See previous section.

Go to 'source/'.

If you use gfortran compiler (default one), just type `make`, the executable 'geoclim.exe' will be created.
If you use a different compiler, you must specify it by typing `make FC=your_compiler`. Note that the configuration is
only made for gfortran, ifort and pgfortran. If you use a different one, the compilation options will not be defined.
See next paragraph for more information.

Note also that if you use the command `build_GEOCLIM`, it will tell you "make command" that was used (make + passed variables).
This command will also be saved in the file "GEOCLIM_environment". It may be useful if you want to recompile the code.

###### Details of Makefile:
The "standard" options to customize the Makefile are:

`make [FC=...] [MODE=...] [ncpath=...] [execut=...] [main_flags="..."] [FFLAGS="..."]`

All those arguments are optional (hence the []).

* FC=...: Fortran compiler. The Makefile is configured for 'gfortran', 'ifort' or 'pgfortran'
  Default: gfortran
* MODE=...: sets the compilation options. 3 options are accepted
    * 'standard': (default), standard check options.
    * 'debug': extra debugging options
    * 'optim': with optimization flags (and less debugging options, like traceback)
* ncpath=...: States the path of the directory where the netCDF-Fortran library is installed. That directory must contain
  the sub-repertories lib/ and include/
  By default, it is '/usr', but you may have installed your netCDF library elsewhere (/usr/local, /usr/local/netcdf, ...)
  It can be inquired with `nc-config --prefix`. The command `build_GEOCLIM` first tries this to get it. 
* execut=...: sets the name of the created executable file. Default is 'geoclim.exe'
* main_flags="...": Override the main compilation flags (all but "-I$ncpath/include -L$ncpath/lib -lnetcdf -lnetcdff").
  The variable 'MODE' becomes useless. Useful if you use a different compiler whose options are not configured in the Makefile.
* FFLAGS="...": Override all compilation flags. The variables 'MODE' and 'ncpath' become useless. Useful if you use a
  different compiler whose options are not configured in the Makefile and who does not support "-I", "-L", "-lnetcdf" or
  "-lnetcdff" options.

In any case, you can check the compilation command that will be used by doing:

`make echo [all the options you want]`

#### Without Makefile
`make` is not necessary to compile the code, but ensure that:
* The pre-compilation configuration is done (see previous section)
* The file 'source/path.inc' exists and contains the line `character(len=*), parameter:: geoclim_path = "..."`
  Where '...' is the path of the GEOCLIM root directory (this file is edited automatically by the Makefile)
* Your compilation command uses the netCDF options. Usually, it must have the options `-I/usr/include -L/usr/lib -lnetcdf -lnetcdff`
  (/usr/lib and /usr/include are the 2 directories where the netCDF library is commonly installed, but it may be elsewhere).
* You use the 'gnu' or 'Fortran 2003' standard for *all* source files (usually `-free` or `-ffree-from`, possibly `-std=gnu` or `-std=f2003`)
* To make sure that the executable is up to date with the source files, do `rm -f *.o *.mod *__gen_mod*` before compiling (or `make clean`)

Goddéris' tip:
* add to your .bash_profile: `alias nebula="ifort -132 -O3 -free -I/usr/local/include/ -L/usr/local/lib/ -lnetcdff -lnetcdf -o"`
* restart your terminal
* compile with `nebula geoclim.exe *.f *.f90`
* compile 4 times (there are 4 levels of nested subroutine calls)
* it is normal to get many errors on the first 3 times this is executed


## After-compilation configuration

#### User interface
The GEOCLIM input/output interface is managed with 2 files (potentially 3):
* 'config/IO_CONDITIONS':
  Main IO interface, provides the paths of all the required files (like the initial conditions files, climate files...), and the
  name of the output files and which variables to output (under which name).
* 'config/cond_p20.dat':
  States the physical and numerical parameters: solver time steps, output writing frequency, duration of run, when to generate restarts,
  biodiversity parameters, acceleration parameters, and volcanic degassing (CO2, SO4, Trapp setting...)
* You can specify in config/IO_CONDITIONS to read the climatic inputs from GCM output files. In that case, the names of those files
  and their variables must be provided in 'config/GCM_input_conditions'.

NOTE: you can find some examples of configuration files in config/templates/


#### Input files
The name of all input files needed by GEOCLIM must be provided in config/IO_CONDITIONS
* Oceanic configuration (basin volumes, surfaces, water circulation between them). Reference files are in INPUT/COMBINE/ref/
* Continental climate: total grid area, land area, temperature (at each CO2 levels) and runoff (at each CO2 levels).
  See examples in INPUT/ERA5/ (ascii files. Area in 1e6 km2, temperature in °C, runoff in cm/y).
  Alternatively, those 4 inputs can be read from GCM outputs (land annual climatology), in netCDF format.
  In that case, the GCM output files info must be provided in config/GCM_input_conditions. The model expect 1 netCDF
  file per CO2 level. It recognizes several units and completes automatic conversions. If the units are not recognized, the default
  units will be assumed, but the user will be asked interactively to validate it. New units and conversion factors can easily be added
  in source/physical_units.f90
* Lithology mask: a netCDF file containing the area fraction covered by each lithology class in every land pixels.
  Alternatively, the user can specify directly in config/IO_CONDITIONS a geographically uniform lithology fraction for each class.
* topography slope: a netCDF file, only read if coupled with DynSoil.
* A special file is also required ("All" file) containing mean climate variables (like global temperature) at each CO2 level.
  See for instance 'INPUT/ERA5/All_ERA5.dat'


#### Initialization files
GEOCLIM also needs an initial condition to start.
* COMBINE: a COMBINE restart file (like restart/geoclim/output.ref), stating the values of the main ocean chemistry variables.
  The name of that file is stated in config/IO_CONDITIONS
* DynSoil (if coupled with): you can use the options 'startup:null' 'startup:eq' (in config/IO_CONDITIONS), to automatically
  generate initial conditions with null regolith (1st case) or regolith at equilibrium with initial climate (2nd case).
  Alternatively, the option 'restart' tells the code to use a DynSoil restart file (defined below).
* Biodiversity module (if coupled with): needs 2 restart file, one for species and the other for chemical variables.
  The name of those 2 files is stated in config/IO_CONDITIONS


## Inputs error handling
GEOCLIM performs tests on the input files before the "main" execution. They are 4 types:
1. Axis mismatch between the input files (for instance, shifted longitude)
2. missing values on continental pixels (continental pixels are defined by the "land area" input variable)
3. invalid value for runoff (negative) and slope (negative or null)
4. sum of all lithology classes differs from 1

By default, the executable interactively asks the user what to do when an error is encountered. This can be problematic when run
as a batch process on a cluster (with no interactive interface). It is possible to pass 4 arguments to the executable, as follows:

`./executable i1 i2 i3 i4`

where i1...i4 are integer numbers, between -1 and 3, and correspond respectively to the 4 kind of errors above-mentioned.
* -1 means 'ask the user interactively' (default)
*  0 means 'abort the execution'
*  1 means 'remove the problematic pixels' (not possible for axis mismatch)
*  2 means 'ignore the issue and continue execution without any change'
*  3 means 'replace the invalid value' (only possible for runoff and slope)

For instance, I recommend `./geoclim.exe 0 1 3 0`, or if you are sure of your axis and lithology mask `./geoclim.exe 2 1 3 2`

NOTE: these tests cannot be done with ascii input files (except the negative runoff test), use that format at your own risk!

Also, errors exist that the code cannot handle, for instance, if one of the input file specified in config/IO_CONDITIONS
does not exist, or if the shape of the netCDF variables does not match the one specified in the code (source/shape.inc).
These will cause the run to crash. The code can, however, handle transposed 2D (x-y) variables in GCM input files, as well as
degenerated (size-1) extra dimensions.


## Executing
`build_GEOCLIM` put the executable file in "executable/". The other methods let it in "source/". It can be run from any
directory, as all the paths in the code are absolute.


## Output
* 2 systematic netCDF files: geoclim and geographic
* 1 netCDF file for DynSoil outputs (if coupled with)
* 1 netCDF file for biodiversity outputs (if coupled with)
* It is possible to automatically convert those outputs in ascii format,
  that option must be entered at the last line of config/cond_p20.dat

You can specify in "config/IO_CONDITIONS" which variable you want to output or not (among the list of "outputable" variables).
To not output a variable, write "#" in place of the output file name.
You can also change the name of the output files, each variable can be placed in a separated file.

The frequency at which outputs are written is specified in "config/cond_p20.dat". There are 3 frequencies: one for COMBINE and
biodiversity output, one for geographic outputs (i.e., continental variables) and one for specific DynSoil outputs.


## Restart
Restarts are created in the output directory. It is a good habit to move them to "restart/.../"

The time for restart generation is specified in "config/cond_p20.dat". By default, this time is set to "-1.", which
is interpreted as "when the run is complete". Note that the restarts are generated 1 time only.


## Killing a run
Sometimes one may want to end a run and create restart files precociously. For instance, if a run has reached the steady-state
sooner than expected and one wants to launch "perturbation runs" from that steady-state.

To do so, simply write the name of the run (as specified at the first uncommented line of "config/IO_CONDITIONS") in the file
"deathnote.txt". It will cause the run to stop and generate restart files (if they were not already created).
You can put as many run names as you want in the deathnote, one by line. Don't forget to erase the names afterward!


## Special runs

#### Equilibrium run
An "equilibrium run" is a run whose transient evolution is of no interest because one only wants to get the geochemical
steady-state (for instance, to start perturbation from that steady-state).
In that case, a couple of things can be modified to shorten the time needed to reach the steady-state, without modifying it.

###### Before compilation: 
Only if you are using DynSoil module in its dynamic version, should you decrease the value of 'scaling_factor' in
"source/dynsoil_physical_parameters.inc". The scaling factor controls the inertia of the regolith, and does not affect
its steady-state. 1 is for a normal regolith. In some places, regolith can take millions of years to reach its steady-state.
To shorten that time, set it to 1d-3. You can put a value as close to zero as you want, it will not generate
any numerical instability. However, it will become useless if the evolution time-scale of the regolith is lower than the
model time step.

Do not forget to put the 'scaling_factor' back at 1 after the run is complete!

Alternatively, you can use the "steady-state" version of Dynsoil. The code will directly compute the analytical steady-state.
This has a lower computation cost, but it has no visible effect, if the asynchronous time step of continental weathering is
high enough.
Note that there will be a slight difference between the analytical steady-state and the numerical one (reached with the
"dynamic" version of DynSoil) simply because of the vertical discretization.

###### After compilation
* Oxygen cycle acceleration: O2 is the species with the longest residence time (~10 Myr). An acceleration coefficient
can be tuned in 'config/cond_p20.dat'. Setting it to 100 is enough to bring it down to Carbon residence time. An excessively
value will cause the model to crash.
* Asynchronous coupling with continental weathering: The standard time-step for continental weathering is 25 years, and
the model spends a significant amount of time on the continental computation, especially at high resolution (1° or less)
and when coupled to DynSoil.
Increasing that time step to 250 years, or 1000 years will hasten the run, only degrading the quality of the transient
evolution. If it is too high, however, it can increase the model time needed to reach the steady-state.

#### Calibration run
To be fully consistent, the model should be recalibrated for each new set of boundary conditions.
The current calibration is done with ERA5 reanalysis fields for temperature and runoff, SRTM slope, and
Hartmann et al. 2013 lithology mask, all at a resolution of 0.5°.
The climate fields of a General Circulation Model will inevitably differ from ERA5 fields, and differ from one GCM to
another, yielding many unique geochemical steady-states. The spatial resolution may also affect the steady-state.

Here are suggested steps to properly recalibrate the model:
* Run a Pre-Industrial simulation (1850 boundary conditions, 1xCO2) with the GCM you intend to use, preferably with
the same set of components and resolution. Retrieve the equilibrium annual climatology.
* Configure GEOCLIM at the given geographic resolution and with 1 CO2 level. It will set the model in fixed CO2 mode.
* For the lithology mask, it is preferable to keep it the same way than you intend to use it for the paleo runs
(uniform or spatially-resolved, with same number of classes), while being consistent with present-day lithology.
* Do a first short run with the pre-industrial forcings to retrieve the silicate weathering flux. If you used DynSoil
in "dynamic" mode, run the model from "startup:eq" with acceleration tuning (see previous paragraph) during 1000-10000 yr.
If you use another set of components, only 1 model time step is needed.
* Get the "total silicate weathering" flux from the outputs and use that value as CO2 degassing flux. It should be 2-6 Tmol/yr.
As today's degassing flux is not well constrained, it is better to tuned it and keep the weathering parameters unchanged.
Note that the degassing flux (specified in "config/cond_p20.dat") is split in 2: Volcanic (continental) and MOR (oceanic).
Doing so, the equilibrium CO2 with Pre-Industrial boundary conditions will be 1 PAL.
* *Phosphorus weathering?*
* With that degassing flux, re-do as many runs as needed to get 1 PAL of atmospheric O2 (use oxygen acceleration coefficient!)
There is no other way than to manually run the model, adjust the parameters if O2 is too low or too high, re-compile, re-run,
and so on, because the organic carbon burial flux also depends on O2 concentration.
I recommend tuning the value of the parameter 'OC_in_rocks' that corresponds to "siliclate sediments", because it is the
most poorly constrained parameter. This parameter (in "source/combine_foam.inc") specifies the mass fraction of petrogenic
carbon in each rock type. A higher value will result in higher kerogen weathering, and thus less oxygen (and vice versa).
A standard value for siliclate sediment is ~1%, though is highly depends on the type of sediment.
Note that the total kerogen weathering carbon flux (= total organic carbon burial) should be ~5 Gmol/yr.
* Though it is of less importance, you can adjust the carbonate weathering constant in "source/cont_weath.f" to have a total
carbonate weathering flux of 13 Gmol/yr. A simple cross-multiplication is sufficient to get the right flux, no need to do
back and forth runs.
However, make sure to do it **before** oxygen tuning, as it may affect the efficiency of organic carbon burial.


## How to generate boundary conditions
...


## Basic code modifications
...


## Visualization?
...


## How to install netCDF-Fortran library

#### Mac OS
Note: the following instructions worked in June 2018, on Mac OS High Sierra 10.13.5

Check the following before completing the steps outlined below:
* check if gfortran is installed by typing `gfortran` onto the command line
    * there may be conflict issues between ifort and gfortran (these instructions are meant for ifort compiler)
    * if the command is recognized, then gfortran is installed, and may need to be removed
* XCode is (probably) required - make sure it is installed, up to date, and the licence has been accepted
    * to accept the licence: `sudo xcodebuild -license`
* if something went wrong with the installation, make sure that the source directories (zlib, hdf5, netcdf) are "fresh"
    * don't use a source directory that has been used before
    * instead, delete the old source directory, redownload and unzip, and attempt installation again
* restarting the computer after step 1 and step 4 may be helpful

1. netCDF-C
    * download zlib:
        * http://www.zlib.net
    * download hdf5
        * https://www.hdfgroup.org/downloads/hdf5/source-code/
    * download netCDF
        * https://github.com/Unidata/netcdf-c/releases/v4.6.1
    * unzip and `cd` into the zlib directory
        * `export ZDIR="/usr/local"`
        * `./configure --prefix=${ZDIR}`
        * `make check`
        * `make install`
    * unzip and `cd` into the hdf5 directory
        * `export H5DIR="/usr/local"`
        * `./configure --with-zlib=${ZDIR} --prefix=${H5DIR} --enable-hl`
        * `make check`
        * `make install`
    * unzip and `cd` into the netCDF directory
        * `export CPPFLAGS="-I${H5DIR}/include"`
        * `export LDFLAGS="-L${H5DIR}/lib"`
        * `./configure`
        * `make check`
        * `make install`
2. netCDF-F
    * download netCDF-F:
        * https://www.unidata.ucar.edu/downloads/netcdf/index.jsp
    * unzip and `cd` into the netCDF directory
        * `export NCDIR="/usr/local"`
        * `export NFDIR="/usr/local"`
        * `export CPPFLAGS="-I${NCDIR}/include"`
        * `export LDFLAGS="-L${NCDIR}/lib"`
        * `./configure --prefix=${NFDIR}`
        * `make check`
        * `make install`

#### With Linux OS
`apt-get install libnetcdff` :)


## Non-required software instructions for installation (Mac OS)
* ifort
    * download the student version
        * https://software.intel.com/en-us/qualify-for-free-software/student
    * click `macOS (Fortran)`
    * follow sign up
    * download and install
        * `source /opt/intel/bin/compilervars.sh intel64`
    * documentation
        * file:///opt/intel/documentation_2018/en/ps2018/getstart_comp_mf.htm
* pyFerret
    * download for python 3.6
        * https://github.com/NOAA-PMEL/PyFerret/releases
    * follow instructions here:
        * https://github.com/NOAA-PMEL/PyFerret#installation-from-prebuilt-targz-file
    * every time a new terminal is opened, you have to run `. ferret_paths.sh` (note the space)
        * recommend adding to `.bash_profile`:
            * `alias enable_pyferret=‘. /path/to/dir/pyferret-7.4-MacOSX-Python-3.6/ferret_paths.sh’`


## Frequent issues

### Errors during compilation

#### incompatibility with previously compiled files
It happens when some files were previously compiled with another compiler, or other options, and are not compatible
with your new compilation command. Type `make clc` in the "source/" directory (or use option `--reset` in `build_GEOCLIM`)
and see if the error still persists.

#### netCDF library
This is by far the most frequent source of error, and they can be hard to detect and solve. They are basically 2 possibilities:

###### library not found
To check if the library is found, try and compile only the file "netcdf_io_module.f90" (that does not use any other source file),
for instance: `make netcdf_io_module.o [your potential Make options]`
If the library is not found, the error message should say that the netCDF module is not found, without even trying to use the
functions. It should look like:

> netcdf_io_module.f90:8:5:
>
>   use netcdf  
>      1  
> Fatal Error: Can't open module file ‘netcdf.mod’ for reading at (1): No such file or directory

This error occurs when no netCDF library exists in the specified path. Check that your compilation command
contains `-I.../include -l.../lib` (for instance, try `make echo`). If it does contain it, the paths `.../include'
and '.../lib' probably do not have netCDF library. Try to find where the library is installed.
This information can normally be obtained with `nc-config --prefix` (note that it is what `build_GEOCLIM` uses).

###### library not recognized
This error is harder to detect. The modules and objects are generally successfully created. The error comes while creating
the executable. The compiler often returns plenty of error messages (which does not make the task any easier). Most of them
would that look like:

> /usr/bin/ld: /tmp/cc7LF9wU.o: in function \`\_\_netcdf_io_module_MOD_put_att_int':  
> netcdf_io_module.f90:(.text+0x8148): undefined reference to \`\_\_netcdf_MOD_nf90_put_att_one_fourbyteint'

This indicates that there is a netCDF library, but the compiler didn't manage to use it. 
It can happen for several reasons:
* Some compilation options are missing. Compilers generally need a specific option to use netCDF, `-lnetcdff` (gfortran),
`-lnetcdf` (ifort). The Makefile put those 2 options among the compilation flags (unless you override the flags with `FFLAGS=...`).
Make sure that you have those 2 options in your compilation command.
Sometimes (with ifort notably) those options need to be *at the very end* of the compilation command, which is what the Makefile
normally does.
* Incompatible library: another possibility is that the netCDF library is not compatible with the compiler. It may be a version
issue, or the fact that it is installed for the wrong compiler (try `nc-config --fc` to check the compiler the library is
configured for). In that case, there is no better solution than reinstalling the netCDF library (or using another Fortran compiler).

#### Fortran fixed-format interpretation
Normally, this type of error should not happen with the Makefile. With plain compilation command, however, one must be careful to
specify free format interpretation (`-ffree-form` with gfortran, `-free -132` with ifort...), as by default. Some compilers consider
files with extension '.f' (or all files) as Fortran 77 format.
This ".f" file interpretation is also an implicit rule in `make`, but this rule is overridden in the current Makefile.

### Error during execution

#### With netCDF input format
i.e., input mode = 'GCM'.

With that format, the code is able to do many compliance checks and detect most error sources (and notify the user).
The code does not check invalid value for area (or land area) and temperature.
Be careful for instance if you define land fraction as a difference of variables, not to generate negative values.

The netCDF-Fortran library may not be able to read file in the most recent netCDF versions, like netCDF4.
The code is meant to read netCDF3 "classic" format input files.
This should not be an issue for GCM outputs, but be careful if you export data (like slope and lithology) in netCDF
format. In python with netCDF4 package, specify "format=NETCDF3_CLASSIC".

#### With ascii input format
i.e., input mode = 'ascii'

The code does not perform any checks besides negative runoff and slope. There can be several sources of errors.
The run will not necessarily crash, it sometimes continues with NaN or Infinity values. If that happens, recompile the code
with debug options (use `MODE=debug` with the Makefile, or `--mode debug --reset` with `build_GEOCLIM`) and re-run it.
It will tell you where the first error happened.

The "standard" GEOCLIM ascii format for geographic fields is:
* values unravelled with increasing latitude and longitude, longitude being the most rapidly varying axis.
* for total area and continental area, one pixel value per line, or values separated by commas or space (works as well)
* for climatic fields (temperature and runoff), **the first value must be the current CO2 level**, then, all the pixels
  values of the field (similarly unravelled), then the next CO2 level, and so on.
  Usually, each line is for one CO2 level, and the values are separated by comma or blank, but it works just as well with
  line breaks, or even one value per line, as long as the order is respected. 
  CO2 levels must be in **decreasing order**.

Slope and lithology mask must be in netCDF format. Note however that you can specify a uniform lithology directly in
"config/IO_CONDITIONS", and that slope file is only needed if DynSoil module is activated.

###### Missing values
It is possible that on some continental points (i.e., points with area > 0), climatic fields have missing value (runoff notably).
Note that the code will always check if there are points with negative runoff.

If the missing value is far enough from the valid range (like -9d33, or 1d36), you should get an error message like:

> Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.

In one of those files:

> 0x55a5bfc51603 in eqcte\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5.2/source/eqcte.f:16

> 0x55fc6d79cb98 in bio\_frac\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5.2/source/bio_frac.f:8

>  0x5578b564025d in carbo\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5.2/source/carbo.f:9

> 0x55c3a02db371 in carbo\_  
>	at /home/piermafrost/GitHub/GEOCLIM5.2/source/carbo.f:24

> 0x55a479b87442 in ocean\_atm\_flu\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5.2/source/ocean_atm_flu.f:15

If the missing value is close from the valid range, there may be no error, simply wrong continental fluxes.

###### Shifted grid
If the formatting of climatic or slope fields is different than the one of area (e.g., flipped latitude axis, or longitude
starting at -180° instead of 0°), the code will likely read missing values, and the same errors than previous paragraph
will happen. However, it may read regular temperature values out of continents, or null runoff, depending on how the input
ascii file handle non-continental points. In that case, you will simply have wrong continental weathering fluxes, which may
be difficult to identify.

###### Units
Ascii files carry no information on variable unit, so the code cannot check it.  
The unit assumed by the code are:
* area: m2
* temperature: °C
* runoff: cm/yr

If the temperature is in Kelvin, you should receive an error message like:

> Program received signal SIGFPE: Floating-point exception - erroneous arithmetic operation.
>
> 0x55a479b87442 in ocean\_atm\_flu\_  
> 	at /home/piermafrost/GitHub/GEOCLIM5.2/source/ocean_atm_flu.f:15

Wrong runoff or area units will generally not trigger a crash, but will generate aberrant continental fluxes:
too high by a factor ~10 (runoff in mm/yr), or too low by a factor that is often 1/10, 1/100, 1d-6, 1d-12.

The output variable "discharge" (water discharge) is a good indicator of wrong runoff or area units, as its order
of magnitude is normally ~4d13 m3/yr.
You can also simply check the climatic and area variables in the geographic output file.

#### Error with COMBINE input data, or initial state
Combine input data (size of oceanic basins, seawater temperature...) is expected to be different for each paleo configuration.
It is possible that a new input dataset may make the model crash, because of error in its generation, or because it is not
compatible with the initial condition.
A reason for that is that the restart file (COMBINE initial condition) gives the model the absolute amount of chemical
species in each basin (in moles). If the volume of one basin is significantly changed from the configuration that the restart
came from, that may generate aberrant concentration of chemical species.

To see if the errors come from Combine input data, try and re-run the model with the reference dataset (files in
"INPUT/COMBINE/ref/").  
If the COMBINE initial condition is solely responsible for the crash, reducing the time step just for the time to the
ocean mixing to dissipate aberrant concentrations (100-1000 years) may solve the problem. The model could then be
run normally from the new restart.


## Notes
Reference for ERA5 climate dataset:
Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I.,
Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2019): ERA5 monthly averaged data on single" levels from 1979 to present.
Copernicus Climate Change Service (C3S). Climate Data Store (CDS). (accessed on 19 Feb 2020)
https://doi.org/10.24381/cds.f17050d7
distributed under Copernicus Products license: https://cds.climate.copernicus.eu/api/v2/terms/static/licence-to-use-copernicus-products.pdf


## Contact
pierre.maffre@normalesup.org


godderis was here
