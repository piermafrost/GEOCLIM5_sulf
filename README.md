# GEOCLIM5.2 - sulf



## Updates from GEOCLIM5:
* Possibility to read area, temperature and runoff inputs directly from GCM annual climatology files,
  with some automatic unit conversions.
* Possibility to automatically create DynSoil initial conditions (null regolith, or at equilibrium with initial pCO2).
* Can specify uniform lithology, without input netCDF file
* Add checks for axis matching, missing points, invalid runoff and slope, invalid lithology, and physical units.
  Interactively ask the users if error detected (by default).
* Reorganize the main IO file ('config/IO_CONDITION'), can now add commented (#) and blank lines
* Standardization of the oceanic temperature input (2 options, "parametric" or "raw ascii input").
* Minor code improvements (eg, Runge-Kutta 4 scheme, biological 13C fractionation formula...)
* Logarithmic or linear interpolation available for CO2 dimension (default: log)
* COMBINE (oceanic) restart in concentration instead of molar amount.
* Implementation of a simplified sulfur cycle (oceanic SO4^2-, sulfide weathering and sulfate-reduction).



## GEOCLIM model in a nutshell
GEOCLIM is cluster of models, more or less adaptable, computing geochemical cycles of several species (C, O...) at geological
timescale.
The core of the model is an ocean-atmosphere chemistry model (advection-reaction) COMBINE. Ocean an atmosphere are discretized
in 10 reservoirs (boxes).
This core is closely associated to an early diagenesis module computing the "output" fluxes (burial of elements in marine sediments)
for each box.
It is also associated to an continental weathering module computing the "input" fluxes. This weathering module is spatially-resolved
(using a geographic mesh grid), its resolution is adaptable, and several options exist for the silicate weathering part.
This triplet is indirectly coupled to a climate model (GCM).
Climate simulations must be run before using GEOCLIM, for a range of CO2 levels. Any climate model can be used, as long as it computes
surface air temperature and continental runoff. Oceanic temperature from the climate model can be used, or parameterized if not
available. Climate fields are then interpolated on the "CO2 dimension", at the current atmospheric CO2 computed by GEOCLIM, and
used to compute continental weathering, and oceanic boxes temperature.
The resolution (ie, the spatial grid) of the continental weathering module must be the same than the GCM.

##### How to run the model:
After downloading the present repository, type `./make_test testname` (testname being one of "ERA5", "CESM", "paleo" and
"ascii"). This command will compile and execute a short GEOCLIM run and compare the output to a reference template.
This allows to verify that the compilation and execution of the model are performed without error, and yield the same
results than reference runs. If not, the command should tell what type of error was encountered (see also section *Frequent issues*
at the end of this file). You could try the 4 tests to make sure everything works as excepted.

If the tests are conclusive, follow those step to create your run:
* Compile the code with `build_GEOCLIM` (specifying the model set of components and resolution). Try `./build_GEOCLIM --help`
for more information. This command create an executable file in 'executable/'
* Configure the run by editing the files 'config/IO_CONDITIONS' and 'config/cond_p20.dat' (name of the run, initial condition,
forcing fields, solver parameters...). This post-compilation configuration must be consistent with the specified set of components
and resolution.
* Run the model with `executable/exec_name` (name of the executable created at the compilation step). Alternatively, and if you
are running the model on a cluster, you can use the files in 'job/' to submit the run as a batch process. 'job/submit_chain.sh'
is a script designed for submitting a series of runs, each new one starting from the end of the previous one. See also section
*Multiple runs and job submission* in *Run GEOCLIM*.



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
The script `make_test` does the complete configuration (pre- and post-compilation), compiles and runs the desired template
(the pre-existing configuration and source files will be restored when the run is complete).
The expected GEOCLIM outputs for those templates are available in "OUTPUT/templates/", so you can check that you obtain
the same results.
Type `./make_test` for more information.



## Useful commands

### `configure.sh`
A script to quickly set up predefined GEOCLIM configurations (pre- and post-compilation), like "ref" (ERA5) or "GFDL".
It replaces some current GEOCLIM files by predefined ones, stored in local template/ directory.
The replaced files are the 3 config files (IO_CONDITIONS, cond_p20.dat and GCM_input_conditions, in config/),
and some Fortran source files, (constante.f90, shape.inc and coupler.inc, in source/).
GEOCLIM is then configured for the given forcings (climate fields, oceanic geometry...) with the corresponding
model parameters calibrated for those forcings.

The "new" configuration files can still be modified after invoking that command, for instance, to run the model
in a paleo configuration using the same GCM and parameterization as a predefined one.

Because this script updates the source code, it should be invoked before compilation. Furthermore, as it performs
the pre-compilation, it is not needed to use the command `build_GEOCLIM` to compile the code, the Makefile (in source/)
is enough. `build_GEOCLIM` can still be used, but the user will have to re-specified the model resolution and set of
components.

### `build_GEOCLIM`
`build_GEOCLIM` is a bash script, meant to edit the needed source files for pre-compilation configuration (data shape,
activated modules, ...) and compile the code. The post-compilation configuration is left to do (which file to use as
input data, names of output and solver parameters, see section "After-compilation configuration").
It uses the Makefile to compile the code.

This command supports a large number of options (specify compiler and options to use...).
Type `./build_GEOCLIM --help` to get detailed information on how to use it.



## Pre-compilation configuration
If you do not use the `build_GEOCLIM` nor `configure.sh` commands, here is the list of configuration steps you should
follow *before compiling the code*, in order to define your run:
* Which modules to use (file 'source/coupler.inc')
* Number of CO2 levels: 'nclimber' parameter in 'source/shape.inc'
* Geographic resolution (must be consistent with input data): 'nlon' and 'nlat' parameters in 'source/shape.inc'
* Number of lithological classes: variable 'nlitho' in 'source/shape.inc'. The model is parameterized for 6 classes, and
  the subroutine 'source/cont_weath.f' explicitly assumes that carbonate rocks is the last lithology class (regardless of
  the number of classes).
  If you want to define more (or less) lithology classes, you should:
    * Leave carbonates as last lithology
    * Update the lithology-dependent parameters in 'source/constante.f90'
    * Update the lithology-dependent parameters in 'source/dynsoil_physical_parameters.f90'
    * Specify which lithology corresponds to "basalt" (variable 'BASALT_LITHO_NUM' in 'source/constante.f90)

  Note that even though 'CaMg_rock' and DynSoil parameters are *not used* if GEOCLIM is not coupled with DynSoil,
  they must be defined consistently with 'nlitho' or the compilation will fail. You may simply fill them with zeros.
* DynSoil vertical resolution: parameter 'nDSlev' in 'source/shape.inc' (must be defined even if not coupled with DynSoil).
  Traditionally 10 vertical levels are used, but it can be changed just by modifying that number.
  If not coupled with DynSoil, its value is not used.
* The number of oceanic basins, as well as the number of geochemical variables, are not meant to be customized, and requires
modification of the source code to be changed. See section *Advanced customization*.
* If you don't compile the code with the Makefile, you should also check that the file 'source/path.inc' exists and states
  the correct path of GEOCLIM root directory (with the line `character(len=*), parameter:: geoclim_path = "..."`). If you do
  use the Makefile, it will be updated automatically.
* A certain number of configuration options are meant for a ecological network module. These options are left for further
development, the ecological module is not available in the present distribution.

All those parameters must be consistent with the initialization and forcing files (initial ocean chemistry, climate files...)



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
`make` is not necessary to compile the code. If you do not use it, make sure that:
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
  acceleration parameters, and volcanic degassing (CO2, SO4, Trapp setting...)
* You can specify in config/IO_CONDITIONS to read the climatic inputs from GCM output files. In that case, the names of those files
  and their variables must be provided in 'config/GCM_input_conditions'.

NOTE: you can find some examples of configuration files in config/templates/


#### Input files
The name of all input files needed by GEOCLIM must be provided in config/IO_CONDITIONS
* Oceanic configuration (basin volumes, surfaces, water circulation between them). Reference files are in INPUT/COMBINE/ref/
* Oceanic temperature: 2 options (specified in config/IO_CONDITIONS). "parametric", the temperature of oceanic boxes are
computed with a parametric function of pCO2. "file_name", the boxes temperature are read in an input ascii file stating the
temperature of each box at each CO2 level. Each row of that file represents a CO2 level, each column states the temperature
of the n-th box, except the first column that states the CO2 level (in ppmv).
* Continental climate: total grid area, land area, temperature (at each CO2 levels) and runoff (at each CO2 levels).
  See examples in INPUT/ERA5/ (ascii files. Area in 1e6 km2, temperature in °C, runoff in cm/y, if INPUT_MODE is 'ascii').
  Alternatively, those 4 inputs can be read from GCM outputs (land annual climatology), in netCDF format.
  Set INPUT_MODE as 'GCM' in config/IO_CONDITIONS. The GCM output files info must be provided in config/GCM_input_conditions.
  The model expect 1 netCDF file per CO2 level. It recognizes several units and completes automatic conversions. If the units are not
  recognized, the default units will be assumed, but the user will be asked interactively to validate it. New units and conversion
  factors can easily be added in source/physical_units.f90 (see section *defining new physical units* in *Further information*).
* Lithology mask: a netCDF file containing the area fraction covered by each lithology class in every land pixels.
  Alternatively, the user can specify directly in config/IO_CONDITIONS a geographically uniform lithology fraction for each class.
* topography slope: a netCDF file, only read if coupled with DynSoil.
* A special file is also required ("All" file) containing mean climate variables (like global temperature) at each CO2 level.
  See for instance 'INPUT/ERA5/All_ERA5.dat'
* Global Mean Surface Temperature (GMST): this input is optional, and is available in 'GCM' input mode only. You can specify
  in config/GCM_input_conditions the netCDF files and variable for 2m temperature *everywhere* (not only on lands). It can be
  stored in different netCDF files than land temperature and runoff (as runoff is very often in *land* GCM output, whose variables
  are defined only on lands, and not on ocean). However, this temperature field **must** be defined on the same grid than the land
  outputs. If nothing is specified, the GMST will not be computed and output by the model.


#### Initialization files
GEOCLIM also needs an initial condition to start.
* COMBINE: a COMBINE restart file (like restart/geoclim/output.ref), stating the values of the main ocean chemistry variables.
  The name of that file is stated in config/IO_CONDITIONS
* DynSoil (if coupled with): you can use the options 'startup:null' 'startup:eq' (in config/IO_CONDITIONS), to automatically
  generate initial conditions with null regolith (1st case) or regolith at equilibrium with initial climate (2nd case).
  Alternatively, the option 'restart' tells the code to use a DynSoil restart file (defined below).



## Run GEOCLIM

### Executable files
`build_GEOCLIM` put the executable file in "executable/". The other methods let it in "source/". It can be run from any
directory, as all the paths in the code are absolute paths.

### Inputs error handling
GEOCLIM performs tests on the input files before the "main" execution. They are 5 types:

1. axis mismatch between the input files (for instance, shifted longitude)
2. missing values on continental pixels (continental pixels are defined by the "land area" input variable)
3. invalid value for runoff (negative) and slope (negative or null)
4. sum of all lithology classes differs from 1
5. units not recognized in netCDF inputs

By default, the executable interactively asks the user what to do when an error is encountered. This can be problematic when run
as a batch process on a cluster (with no interactive interface). It is possible to pass 5 arguments to the executable, as follows:

`./executable i1 i2 i3 i4 i5`

where i1...i5 are integer numbers, between -1 and 3, and correspond respectively to the 5 kind of errors above-mentioned.
* -1 means 'ask the user interactively' (default)
*  0 means 'abort the execution'
*  1 means 'remove the problematic pixels' (not possible for axis mismatch or units not recognized)
*  2 means 'ignore the issue and continue execution without any change'
*  3 means 'replace the invalid value' (only possible for runoff and slope)

For instance, I recommend `./geoclim.exe 0 1 3 0 0`, or if you are sure of your axis, lithology mask and units `./geoclim.exe 2 1 3 2 2`

NOTE: these tests cannot be done with ascii input files (except the negative runoff test), use that format at your own risk!

Also, errors exist that the code cannot handle, for instance, if one of the input file specified in config/IO_CONDITIONS
does not exist, or if the shape of the netCDF variables does not match the one specified in the code (source/shape.inc).
These will cause the run to crash. The code can, however, handle transposed 2D (x-y) variables in GCM input files, as well as
degenerated (size-1) extra dimensions.

### Output
* 2 systematic netCDF files: geoclim and geographic
* 1 netCDF file for DynSoil outputs (if coupled with)
* It is possible to automatically convert those outputs in ascii format,
  that option must be entered at the last line of config/cond_p20.dat

You can specify in "config/IO_CONDITIONS" which variable you want to output or not (among the list of "outputable" variables).
To not output a variable, write "#" in place of the output file name.
You can also change the name of the output files, each variable can be placed in a separated file.

The frequency at which outputs are written is specified in "config/cond_p20.dat". There are 3 frequencies: one for COMBINE
output, one for geographic outputs (ie, continental variables) and one for specific DynSoil outputs.

### Restart
Restarts are created in the output directory. It is a good habit to move them to "restart/.../". Automatic launching script,
like submit_chain.sh (in job/) will automatically move the restart files in that directory.

The time for restart generation is specified in "config/cond_p20.dat". By default, this time is set to "-1.", which
is interpreted as "when the run is complete". Note that the restarts are generated 1 time only.

### Killing a run
Sometimes one may want to end a run and create restart files precociously. For instance, if a run has reached the steady-state
sooner than expected and one wants to launch "perturbation runs" from that steady-state.

To do so, simply write the name of the run (as specified at the first uncommented line of "config/IO_CONDITIONS") in the file
"deathnote.txt". It will cause the run to stop and generate restart files (if they were not already created).
You can put as many run names as you want in the deathnote, one by line. Don't forget to erase the names afterward!

### Multiple runs and job submission
'submit_chain.sh' (in the repertory job/) is a bash script for automatically launching a series of GEOCLIM runs.
A series of runs are runs that have exactly the same configuration (except for their timesteps, starting and stopping times),
each one starts from the end of the previous one (the very first one starts from the initial condition given by the user).
This is useful for runs with an initial perturbation, requiring a short timestep, but whose long-term evolution (after
the adaptation to the perturbation) can be computed with a longer timestep.
In addition, it offers to possibility to submit the GEOCLIM run as batch processes (jobs), which is required on clusters.
Clusters usually have a time limitation for jobs, which makes the automatic resubmission (series of run) helpful.
Finally, this script provides a security for conflicting access to the configuration files, that is helpful for running
several independent runs in parallel.

Practically, the script 'submit_chain.sh' works in pair with a second script (usually, 'run_geoclim.sh'). The main script
('submit_chain.sh') "submits" the second one (either executes it, or submits it with the cluster submission command), that
actually run the geoclim model, and call the first script back when the run is completed.
The main script does all the configuration, and move the restarts. The second is only for running the GEOCLIM executable,
but must be configured for the current cluster (whereas the main one is a bash script meant to be executed directly).

When using 'submit_chain.sh', the pre-compilation configuration must be done (and the code compiled). If you wants to launch
in parallel several runs that need different pre-compilation configurations, save as many different GEOCLIM executable files.
Here is the list of options that can be customized with 'submit_chain.sh':
* The name of the run (for a series of runs, suffix '_1', '_2'... are automatically added).
* The submission command (cluster-dependent) and the name of the running script (usually, 'run_geoclim.sh').
* The name of the GEOCLIM executable file.
* GEOCLIM (COMBINE) and DynSoil initial condition.
* Stopping (and restarting) times. Note: The "first" starting time is given by config/cond_p20.dat, and should normally be 0.
* The different model timesteps and printing timesteps.
* The job log file.
* The name of GEOCLIM main configuration file (normally, config/IO_CONDITIONS. In case extra configuration customization
is needed. Usually, keeping the default one is sufficient).

The script is designed for parallel runs. It edits the configuration files and ensures there is no conflict.
Once you have submitted one run (series of run), you can safely edit the file 'submit_chain.sh' and submit a second run (series of
runs). The script will tell you if a run is waiting to access the configuration files.
If you need to do configuration modifications not available in 'submit_chain.sh', the safest way is to create a new config file
"IO_CONDITIONS" and to tell 'submit_chain.sh' to use it (note: the name of the other config files, like cond_p20.dat, are stated
in the main one). Remember that if you edit any of the configuration files, it will impact **all** the series of run that have
been launched. When all the runs are completed, the original configuration files will be reinstated.


### Special runs

#### Equilibrium (accelerated) run
An "equilibrium run", or accelerated run, is a run whose transient evolution is of no interest because one only wants to get
the geochemical steady-state (for instance, to start perturbation from that steady-state).
In that case, a couple of things can be modified to shorten the time needed to reach the steady-state, without modifying it.

###### Before compilation: 
Only if you are using DynSoil module in its dynamic version, should you decrease the value of 'scaling_factor' in
"source/dynsoil_physical_parameters.f90". The scaling factor controls the inertia of the regolith, and does not affect
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
* Oxygen cycle acceleration: O2 has a residence time of ~8 Myr, so it requires around 20 Myr to reach equilibrium.
An acceleration coefficient can be tuned in 'config/cond_p20.dat'. Setting it to 100 is enough to bring the equilibration
time down to Carbon residence time. An excessively value will cause the model to crash.
* Sulfur cycle acceleration: Similarly to oxygen, an acceleration coefficient can be tuned in 'config/cond_p20.dat' to reduce
the time needed for sulfur cycle to reach equilibrium (the residence time of sulfur is ~30 Myr). 100 is a good value.
* Asynchronous coupling with continental weathering: The standard time-step for continental weathering is 25 years, and
100 years for DynSoil module (if activated).
The model spends a significant amount of time on the continental computation, especially at high resolution (1° or less)
and when coupled to DynSoil.
Increasing that time step to 250 years, or 1000 years will hasten the run, only degrading the quality of the transient
evolution. If it is too high, however, it can increase the model time needed to reach the steady-state.
Moreover, the steady-state weathering flux of DynSoil module (in its "dynamic" version) are actually dependent of the DynSoil
timestep, because of numerical accuracy. A longer timestep will result in slightly lower weathering flux (ie, higher equilibrium
CO2). For instance, increasing the timestep from 100 year (default) to 10000 years cause the CO2 to rise by a ~10 ppmv.

#### Fixed CO2 run
This is a special case of model configuration. If there is only 1 CO2 level (nclimber=1), it will be run in "fixed CO2 mode".
This means the amount of CO2 in the atmospheric reservoir will be held constant at the value of the unique CO2 level, whatever
the carbon fluxes. The concentration of the various forms of carbon in the other reservoirs will adjust freely to the atmospheric
concentration (ocean-atmosphere diffusion) and to the carbon fluxes.
In other words, in fixed CO2 mode, the mass balance is not respected for carbon.

This is useful for calibration runs where one wants to hold the atmospheric CO2 constant and adjust the degassing to balance the
silicate weathering flux.

#### Run with locked geochemical cycles
This can be useful if one wants to investigate the behavior of inorganic carbon cycle only, without the feedback of the
other cycles, while still respecting mass-balance. 2 geochemical cycles can be "locked": the sulfur cycle, and the oxygen cycle.
The model "lock" a cycle be imposing that the sources balance the sink at each timestep. More specifically, the sinks (oceanic
processes) are computed freely, and the sources (continental weathering) are force to match the sinks.
* For the sulfur cycle, the sulfuric silicate weathering is still compute freely, and the sulfuric carbonate weathering is adjusted
so that the sum of the two matches the sulfate reduction (the release of H2SO4 is set to 0).
* For the oxygen cycle, the kerogen weathering is adjusted so that when added to sulfide weathering, it matches the organic carbon
burial (whether or not the sulfur cycle is locked).

The mass balance is still respected, which means those modified fluxes affect the other geochemical species (carbon, alkalinity...)

To lock one or several cycle, set the value `.true.` of the corresponding parameters in source/coupler.inc (*before compilation*),
or use the options `--lock OS` (O for oxygen, S for sulfur, a single one works as well) in `build_GEOCLIM`.



## Further information

### Calibration procedure
To be fully consistent, the model should be recalibrated for each new set of boundary conditions.
The current calibration is done with ERA5 reanalysis fields for temperature and runoff, SRTM slope, and
Hartmann et al. 2013 lithology mask, all at a resolution of 0.5°. A second calibration is available for the GFDL boundary
conditions.
The climate fields of a General Circulation Model will inevitably differ from ERA5 fields, and differ from one GCM to
another, yielding many unique geochemical steady-states. The spatial resolution may also affect the steady-state.

Here are suggested steps to properly recalibrate the model:
* Run a Pre-Industrial simulation (1850 boundary conditions, 1xCO2) with the GCM you intend to use, preferably with
the same set of components and resolution. Retrieve the equilibrium annual climatology.
* Configure GEOCLIM at the given geographic resolution and with 1 CO2 level. It will set the model in fixed CO2 mode.
If you are using the 'GCM' input mode, you will need to remove (or comment) the lines stating the netCDF inputs that are
not at 1xCO2 in config/GCM_input_conditions. If you are using the 'ascii' input mode, you will need to remove the not-1xCO2
inputs *in the ascii files* (or create new ascii files with only 1xCO2 inputs).
* For the lithology mask, it is preferable to keep it the same way than you intend to use it for the paleo runs
(uniform or spatially-resolved, with same number of classes), while being consistent with present-day lithology.
* Do a first short run with the pre-industrial forcings to retrieve the silicate weathering flux. If you used DynSoil
in "dynamic" mode, run the model from "startup:eq" with acceleration tuning (see previous paragraph) during 1000-10000 yr.
If you use another set of components, only 1 model time step is needed.
* Get the "total silicate weathering" flux from the outputs and use that value as CO2 degassing flux. It should be 2-6 Tmol/yr.
As today's degassing flux is not well constrained, it is better to tuned it and keep the weathering parameters unchanged.
Note that the degassing flux (specified in "config/cond_p20.dat") is split in 2: Volcanic (continental) and MOR (oceanic).
Doing so, the equilibrium CO2 with Pre-Industrial boundary conditions will be 1 PAL.
* Although it is not strictly necessary, you may want to adjust the parameters of Phosphorus and carbonate weathering to get the
desired flux. Phosphorus weathering will impact the oxygen levels. Carbonate weathering has no impact on equilibrium CO2, and
virtually no impact on O2 (though it may affect the biological pump). However, it directly impacts the oceanic DIC, Calcium and
alkalinity. Phosphorus weathering parameters (ie, P amount in source rocks: P_rock, P2C_ker and P2C_carb) are defined in
source/constante.f90. Carbonate weathering should be modified directly in cont_weath.f. In both cases, a simple cross-multiplication
is sufficient to get the right flux.
* Finally, re-do as many runs as needed to get 1 PAL of atmospheric O2 and 29 mol/m3 of mean oceanic sulfate (using acceleration
coefficients will help). There is no other way than to manually run the model, adjust the parameters if O2 and sulfate are too
low or too high, re-compile, re-run, and so on. I recommend tuning the value of the parameter 'OC_in_rocks' (in source/constante.f90)
that corresponds to "silicate sediments", because it is the most poorly constrained parameter. This parameter specifies the mass
fraction of petrogenic carbon in each rock type. A higher value will result in higher kerogen weathering, and thus less oxygen (and vice versa).
A standard value for silicate sediment is ~1%, though is highly depends on the type of sediment. For the sulfur cycle, the parameter
'Sulf_rock' (amount of sulfide in source rocks, still in source/cont_weath.f) is controlled by the S:C ratio, and determine the
sulfide weathering flux. With the acceleration parameters at 100, the model should be run for 2-5 Myr to have an idea of the equilibrium
O2 and SO4^2-
* To be more accurate, the oceanic alkalinity, DIC and the O2 gradient can be checked. If needed, they can be adjusted by modifying the
parameters controlling ocean chemistry, that are defined in source/constante.f90.

### Defining new physical units
When reading netCDF inputs (in "GCM" input mode), the code read the attribute "units" (a string) of the netCDF variables. If the "units"
string matches a defined ones, the corresponding conversion is performed to set the variable into the model's reference units.
It may be needed to add new unit string if the code do not recognize a given netCDF input file (for instance, in my experience, there are
as many runoff units as there are GCM). The string has to match exactly to be recognized (space and case sensitive).

Physical units are defined in the source code in 'source/physical_units.f90'. To define a new one, go to the desired variable (eg,
'temperature_units'), increment the variable 'naccepted' by 1, then go to the last line of 'accept_unit' definition, and add a line defining
'accept_unit(n)%string' (your unit string) and 'accept_unit(n)%conversion' (conversion factor and offset), 'n' being the number of your
newly defined units. The rule for converting variable into the reference unit is 'ref_unit_var = factor\*var + offset'.

### Generate oceanic temperature ascii input file from GCM outputs
This is still experimental, but the Python file 'make_oceanic_input.py' in 'preproc/python/' is meant to read oceanic temperature from GCM
netCDF output for each CO2 level, as well as the ocean grid definition (ie, longitude, latitude and depth), average the temperature on GEOCLIM
(COMBINE) oceanic boxes and generate a ascii input file readable by the model (1 CO2 level per row, CO2 value + boxes temperature on each column).
See description in the Python file.

### Other boundary conditions
##### Note on climatology average
The continental climate (temperature and runoff), as well as the oceanic temperature, that are considered by GEOCLIM model, are supposed to be
annual mean variables, averaged over many years (climatology average. eg, 30 years). The model does not take into account the seasonal cycle, or
other form of climatic variability.

##### Oceanic boundary conditions
Although it is the only GCM and CO2 dependent one, seawater temperature is not the only oceanic boundary condition considered by
the model. The rest of the boundary conditions are defined in 'INPUT/COMBINE/ref/' (the name of the input files can actually be
changed in 'config/IO_CONDITIONS').

They consist of:
* boxes (ie, oceanic basin) definition: a series of ascii files define with 0/1 which boxes are surface (indice_surface.dat), deep
(indice_deep.dat), thermocline (thermocline.dat), polar (indice_polar.dat), epicontinental (indice_epicont.dat), as well as which
receive continental (ie, riverine) fluxes (apport_ct.dat) and which ones are connected to the seafloor (indice_sedi.dat).
* boxes geometry: a series of ascii files define the volume ('oce_vol.dat', in 1e15m3), top area ('oce_surf.dat') and bottom
area ('surf_sedi.dat', both in 1e15m2) of each box. The volume is actually define for each box and each "internal" variable,
because the variables representing isotopic variables are treated differently than the others, and the volume must be 1e-15 (=> 1)
for them. The volume of the atmosphere must also be 1e-15 because the atmospheric variables are in molar amount, and not in
concentration.
* box pressure: the file 'press_box.dat' define the mean pressure of each box (used for chemical constants). Technically, there are
2 values per line, temperature and pressure, but the first value is ignored (as temperature is read in another file).
* particle sinking rate: the files 'fsink.dat' and 'fsink_inorg.dat' define the sinking parameter of organic and inorganic particles
(respectively) on each box. That parameter is technically the fraction of particles that leave the box by sinking per unit of time
(ie, per year).
* oceanic circulation: 'exchange_2.dat'. Unless the other files that are single-column ascii files (one box per line), this one
is a "2D matrix" ascii file. Each point {i,j} (ie, line i, column j) represents the flux of water (in Sv, ie, 1e6m3/s) going from
the box j (column) to the box i (line). Hence, *the sum of each line must yields the same vector than the sum of each column*. 
The advection of the geochemical species in computed as 'water flux time concentration'.

Though the oceanic circulation is CO2-dependent, it is not taken into account in the code (the model assumes this dependence to be
negligible).

### Basic code modification
##### Model parameters
Most of the empirical parameters of GEOCLIM are defined in 'source/constante.f90'. Those parameters concern the oceanic and diagenesis
components (COMBINE module) and the continental weathering parameters, with the exception of DynSoil module parameters.
However, there are still some parameters defined, or written "in hard" into the subroutines.

The chemical equilibrium constants, for oceanic chemistry (like carbonate speciation), are computed dynamically in 'source/eqcte.f90',
because they depends on temperature, pressure and salinity. Those relationships are better constrained by thermodynamics, and less
susceptible to be modified.

All the parameters of DynSoil module are defined in 'source/dynsoil_physical_parameters.f90'

##### Output additional variables
There are a certain number of predefined output variables that the user can choose to output or not simply by editing the main
configuration file 'config/IO_CONDITIONS'.
However, to output a variable that is not in that predefined list, here are the steps to follow.

###### COMBINE variable
ie, oceanic variable, that have a value for each COMBINE box (like salinity), or a single value (for instance, atmospheric
variable, or continental flux).

First, you need to know the name of the variable *in the Fortran source code*, or the way to compute it.

* Increment by 1 the parameter 'nGEOoutvar' in 'source/combine_foam.inc'.
* Increment by 1 the parameter 'nvar' in 'source/geoclim_create_output.f90'. It must have the same value than 'nGEOoutvar' (this
subroutine does not use combine_foam.inc).
* Add a value in the vector parameter 'indice_boxdef' in both 'source/geoclim_create_output.f90' and 'source/geoclim_write_output.f90'
(that vector must be the same in both files). The value must be 0 is the output variable is a scalar (single value), or 1 if it has
a value for each COMBINE box.
* In the main configuration file (config/IO_CONDITIONS), in the section "OUTPUT CONDITIONS", "GEOCLIM", add a line (eg, copy and
paste the last line) stating the name of the variable *in the netCDF output file*, its units, fill-value and description.
* In 'source/geoclim_write_output.f90', at the end of the section "write variables", (just before the section "output file
closing"), add a block of 3 lines (eg, copy and paste the last 3 lines) that looks like:
> !  
> i = 107  
> if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i), (/real(...)/),   (/nt(i)/), (/1/) ) 

With the *code* variable instead of '...'. If that variable is defined on each box, the line should be a little different:
> if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i), real(...(:)), (/1,nt(i)/), (/nbasin,1/)  )

The output variable can be computed directly in that outputting line.

###### Geographic variable
ie, 2D geographic field (for instance, runoff). This is more complex because the variable can also be defined on the lithology
dimension (it is then a 3D field).

* Increment by 1 the parameter 'nGEOGoutvar' in 'source/combine_foam.inc'.
* Increment by 1 the parameter 'nvar' in 'source/geographic_create_output.f90' and 'source/geographic_write_output.f90'.
* If the output variable is defined on the lithology dimension (currently, only the variables # 11 and 12 are lithology-defined):
    * add the variable number to the loop lines 209-225 of 'source/geographic_create_output.f90'.
    * add the variable number to the loop lines 307-125 of 'source/geographic_create_output.f90'.
    * make the loop lines 317-325 of 'source/geographic_create_output.f90' go from 13 to *nvar-1* instead of from 13 to nvar, as
      the last output variable (ie, the one just added) is defined lithology-defined.
* If the output variable is not defined on the lithology dimension, nothing more should be changed in 'source/geographic_create_output.f90'.
* In the main configuration file (config/IO_CONDITIONS), in the section "OUTPUT CONDITIONS", "GEOGRAPHIC", add a line (eg, copy
and paste the last line) stating the name of the variable *in the netCDF output file*, its units, fill-value and description.
* In 'source/geographic_write_output.f90', at the end of the section "write variables", (just before the section "output file
closing"), add a block of 4 lines (eg, copy and paste the last 4 lines) that looks like:
> i = 15  
> if (filenum(i)>0) then  
>   call put_var_real2D( fileid(i) , varid(i) , real(reshape(..., shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )   
> end if

With the *code* variable instead of '...'. If that variable is lithology-defined on each box, the line should be:
>   call put_var_real3D( fileid(i) , varid(i) , real(reshape(..., shape=(/nx,ny,nlit/), order=(/3,1,2/))) , (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )

* If that output variable cannot be computed from the variable present in 'source/geographic_write_output.f90', declare the code
variable you need as input arguments of the subroutine 'geographic_write_output', and add them in 'source/printf.f', where the
subroutine is called.

###### DynSoil variable
Similar to geographic variable, but output only if DynSoil module is activated, and can be defined on lithology and/or vertical
dimension.

* Increment by 1 the parameter 'nDSvar' in 'source/combine_foam.inc'.
* Increment by 1 the parameter 'nvar' in 'source/dynsoil_create_output.f90' and 'source/dynsoil_write_output.f90'.
* Look on 'source/dynsoil_create_output.f90' for all the specific cases for variable defined (or not) on lithology and vertical
dimensions.
* In the main configuration file (config/IO_CONDITIONS), in the section "OUTPUT CONDITIONS", "DYNSOIL", add a line (eg, copy and
paste the last line) stating the name of the variable *in the netCDF output file*, its units, fill-value and description.
* In 'source/dynsoil_write_output.f90', at the end of the section "write variables", (just before the section "output file
closing"), add a block of 4 lines (eg, copy and paste the last 4 lines) similarly to 'geographic_write_output.f90', paying special
attention on which dimension the variable is defined on.
* Similarly to geographic variables, if needed, declare new input argument and add them in 'source/printf.f', where the subroutine
is called.

### Advanced customization
##### Change oceanic basin definition
The definition of COMBINE basins is supposedly entirely controlled by the input files (in 'INPUT/COMBINE/ref/', see section
*Oceanic boundary conditions*), and could be flexible. However, there are several points in the source code that make it tricky
to customize:
* The number of boxes (parameter 'nbasin', standard number is 10) is defined in 'source/shape.inc', and *systematically set as 10*
by the script `build_GEOCLIM` (line 660). Make sure to update those if you want to change the number of boxes.
* Some routines (mostly 'creades.f') explicitly assume that the boxes are ordered **top to bottom**: for every non-surface boxes k,
the fluxes of sinking particles is taken from the box k-1.
* The atmospheric box must be **the last box** (this is how it is identified by the code routines). Moreover, the indices corresponding
to the atmospheric box **must be kept as it is currently** (1 in 'apport_ct.dat', 0 for all the others), or the loops on boxes
will crash.
* Some routines may assume that there is only 1 box that receives the continental fluxes (currently, epicontinental surface box).
If you divide that box in several ones, there is a risk that the continental fluxes will be multiplied.
* Update accordingly all the box definition files, **and the COMBINE restart file**, that contains a list of value for all variables
in all boxes (the "box" dimension being the most rapidly varying one). Pay particular attention to the input file 'oce_val.dat'
that define the volume for all boxes and variables (with still "box" dimension being the most rapidly varying one). Do keep a
volume of 1d-15 for atmosphere box and for all isotopic variables.

##### Add a new geochemical species
Currently, there are 20 main geochemical species, that are species whose oceanic advection is computed, as well as ocean-atmosphere
exchange (only for O2 and CO2) and sinking for particulate variables. They are:
1. DIC (Dissolved Inorganic Carbon)
2. Alkalinity
3. Phosphate (PO4^3-)
4. Calcium (Ca^2+)
5. Strontium (Sr^2+)
6. Sr in PIC (Particulate Inorganic Carbon)
7. POP (Particulate Organic Phosphorus)
8. PIP (Particulate Inorganic Phosphorus)
9. POC (Particulate Organic Carbon)
10. PIC (Particulate Inorganic Carbon)
11. Oxygen (dissolved + atmospheric)
12. atmospheric CO2 (value set to 0 for oceanic boxes)
13. delta 13 C of DIC
14. delta 13 C of PIC
15. delta 13 C of POC
16. delta 13 C of atmospheric CO2 (value set to 0 for oceanic boxes)
17. 87 Sr / 86 Sr isotopic ratio
18. Not attributed
19. Not attributed
20. Sulfate (SO4^2-)

The variables are in concentration (mol/m3) in the oceanic boxes and in amount (mol) in the atmospheric box, except for the isotopic ones.

Here are the instruction to add a new main ocean-atmosphere:

* Update the total number of main variables (parameter 'nvar_real' in 'source/combine_foam.inc')
* Add a section in 'source/creades.f', which is the subroutine defining the input-output fluxes due to chemical reaction, continental
input, sinking (for particulate variables only) and sedimentation on seafloor. This net local I/O rate is the code variable 'R'.
* Add a section in 'source/derivs.f' defining the derivative, that is the sum of 'R' and the oceanic advection **for dissolved
variables only**. Particulate variables are not affected by oceanic advection, and the advection part for isotopic ratio (of
dissolved variables) is including in 'R' in 'source/creades.f', because its mathematical expression is different from concentration
variables.
* Update COMBINE restart file and the input file 'oce_vol.dat' (see previous section *Change oceanic basin definition*).
* Update whatever routine needed to compute the geochemical fluxes of that new variable.

### Visualization
A couple of external scripts are designed for the visualization of GEOCLIM output, in 'visualization/python/' (Python scripts) and
'visualization/jnl/' (Ferret scripts)

'visualization/python/plot_final_state.py' draw the main oceanic profiles and geochemical fluxes at the end of a run and save them
in two pdf files 'final_fluxes--\*.pdf' and 'final_ocean_chemistry--\*.pdf'.
Usage: `python plot_final_state.py geoclim_output_file_path`.



## Technical notes

### Install ifort compiler on Mac OS
Note: the following instructions worked in June 2018, on Mac OS High Sierra 10.13.5

* download the student version
    * https://software.intel.com/en-us/qualify-for-free-software/student
* click `macOS (Fortran)`
* follow sign up
* download and install
    * `source /opt/intel/bin/compilervars.sh intel64`
* documentation
    * file:///opt/intel/documentation_2018/en/ps2018/getstart_comp_mf.htm

### How to install netCDF-Fortran library

#### With ifort compiler, on Mac OS
Check the following before completing the steps outlined below:
* check if gfortran is installed by typing `gfortran` onto the command line
    * there may be conflict issues between ifort and gfortran (these instructions are meant for ifort compiler)
    * if the command is recognized, then gfortran is installed, and may need to be removed
* XCode is (probably) required - make sure it is installed, up to date, and the license has been accepted
    * to accept the license: `sudo xcodebuild -license`
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
The simplest way is with apt-get: `apt-get install libnetcdff`
This will ensure the compatibility with the installed Fortran compiler.

### Non-required software instructions for installation (Mac OS)
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
ie, input mode = 'GCM'.

With that format, the code is able to do many compliance checks and detect most error sources (and notify the user).
The code does not check invalid value for area (or land area) and temperature.
Be careful for instance if you define land fraction as a difference of variables, not to generate negative values.

The netCDF-Fortran library may not be able to read file in the most recent netCDF versions, like netCDF4.
The code is meant to read netCDF3 "classic" format input files.
This should not be an issue for GCM outputs, but be careful if you export data (like slope and lithology) in netCDF
format. In python with netCDF4 package, specify "format=NETCDF3_CLASSIC".

#### With ascii input format
ie, input mode = 'ascii'

The code does not perform any checks besides negative runoff and slope. There can be several sources of errors.
The run will not necessarily crash, it sometimes continues with NaN or Infinity values. If that happens, recompile the code
with debug options (use `MODE=debug` with the Makefile, or `--mode debug --reset` with `build_GEOCLIM`) and re-run it.
It will tell you where the first error happened.

The "standard" GEOCLIM ascii format for geographic fields is:
* values unraveled with increasing latitude and longitude, longitude being the most rapidly varying axis.
* for total area and continental area, one pixel value per line, or values separated by commas or space (works as well)
* for climatic fields (temperature and runoff), **the first value must be the current CO2 level**, then, all the pixels
  values of the field (similarly unraveled), then the next CO2 level, and so on.
  Usually, each line is for one CO2 level, and the values are separated by comma or blank, but it works just as well with
  line breaks, or even one value per line, as long as the order is respected. 
  CO2 levels must be in **decreasing order**.

Slope and lithology mask must be in netCDF format. Note however that you can specify a uniform lithology directly in
"config/IO_CONDITIONS", and that slope file is only needed if DynSoil module is activated.

###### Missing values
It is possible that on some continental points (ie, points with area > 0), climatic fields have missing value (runoff notably).
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
If the formatting of climatic or slope fields is different than the one of area (eg, flipped latitude axis, or longitude
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



## Notes and acknowledgements
Reference for ERA5 climate dataset:
Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I.,
Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2019): ERA5 monthly averaged data on single" levels from 1979 to present.
Copernicus Climate Change Service (C3S). Climate Data Store (CDS). (accessed on 19 Feb 2020)
https://doi.org/10.24381/cds.f17050d7
distributed under Copernicus Products license: https://cds.climate.copernicus.eu/api/v2/terms/static/license-to-use-copernicus-products.pdf



## Contact
pierre.maffre@normalesup.org


godderis was here

