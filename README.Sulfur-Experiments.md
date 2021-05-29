# Sulfide weathering experiments from Maffre et al., submitted to GRL (2021)

Information about GEOCLIM runs for the sulfide weathering perturbation experiments

## Run outputs

#### Where to find the output

The outputs of the sulfide weathering simulations are stored in "OUTPUT/sulfur\_experiments/".
Only the COMBINE outputs (i.e., ocean and atmosphere "box" variables and integrated fluxes) are available.
The geographic outputs require too much memory be be stored in GitHub.

#### Run nomenclature

Several run outputs are stored in the repertory. All their names follow the same structure:  
`geoclim_output.runname_?.nc`  
Where "runname" is the name of the GEOCLIM run, and "?" is integer number, denoting the subdivision of a single run.

Simulations with abrupt perturbation where run in 3 parts: 0--100kyr ("\_1"), 100kyr--1Myr ("\_2") and 1--50Myr ("\_3"),
each run restart from the end of the previous run. The asynchronous timestep of continental weathering and DynSoil modules
and of output printing were increase for each following "sub-run", everything else unchanged.

##### Correspondance between run names and simulation names in the article

The names of perturbation run with increase sulfide (pyrite) weathering start by `.PyrW+50`, the ones with proportional
increase of sulfide and petrogenic carbon (kerogen) weathering startss by `.PyrKerW+10-noP` for run with phosphorus weathering
left unchanged, and `.PyrKerW+10-P`.
In addition, runs where conducted with petrogenic carbon weathering increased by 50% (everything else unchanged), and with *oxygen
cycle, sulfur cycle and DynSoil module accelerated, for fast reaching of steady-state*.
Those run names start by `.ker+50`, and do not have the "\_1" ending.

For simulations with increase sulfide weathering, they are 5 possibilities for the fate of the additionnal generated H2SO4. They
are denoted by the suffixes:
* `-Carb`: "carbonate sulfuric weathering perturbation" (presented in main text)
* `-Sil`: "silicate sulfuric weathering perturbation" (presented in main text)
* `-Carb-trdf`: "carbonate trade-off" (presented in the Supporting Information)
* `-Sil-trdf`: "silicate trade-off" (presented in the Supporting Information)
* `-H2SO4`: "H2SO4 release" (presented in the Supporting Information)

Most of the runs were repeated with different strength of oxygen feedback. The corresponding following suffixes were added to the
run names:
* `_noO2fdbk`: case "no-feedback"
* `_noPfdbk-red`: case "feedback-2"
* `_noPfdbk`: case "feedback-1"
* `_Phyd-fdbk`: case "feedback+1"
* `_Phyd-fdbk-plus`: case "feedback+2"
* `_Phyd-land-fdbk`: case "feedback+3"
* `_Phyd-land-fdbk-plus`: case "feedback+4"

For the simulations where the perturbation is set up progressively over 40Myr, the suffix `-prog` is added between the "perturbation
type" and "oxygen feedback" suffixes. Those simulations were conducted in one run, without adapting the timesteps, so only the "\_1"
output exist.

Additional simulations, presenting in the SI, where run. Their run name suffixes are:
* `_cst-oceT`: run with fixed oceanic temperature
* `_cst-CO2`: run with fixed CO2

Finally, one simulation (carbonate sulfuric weathering) was conducted with accelerated oxygen cycle, sulfur cycle and DynSoil module,
to obtain the multi-cycle steady-sate following the perturbation. This run name has the suffix `-fast`, and no "\_1" ending.

**Example of run name**:  
`geoclim_output.PyrW+50-Sil_Phyd-land-fdbk_2.nc`  
2nd run (100kyr--1Myr) of the simulation with the abrupt perturbation "sulfuric silicate weathering" and enhanced oxygen feedback strength
case "feedback+2".

## How to reproduce the runs

Here are the steps to follow to reproduce the runs presented in the article:

1. **Configure the model for the GFDL forcings**:

    The model is currently configured for ERA5 forcings. To link the paths to the GFDL forcings and use the GFDL-calibration parameters,
    type (in the main repertory):  
    `./configure.sh GFDL` 

2. **Set the perturbation in the source code**:

    Go to the source directory (`cd source/`).
    The code is now calibrated for the GFDL pre-industrial steady-state. Running the code without any modification will leave all
    geochemical species constant at their steady-state values.

    The modifications made for the different "perturbation" simulations are:

    * Sulfide (pyrite) weathering perturbation:

        In the source file "cont_weath.f", divided in 3 steps:
        * lines 312-323: uncomment one of the lines `faddsulfw = ...` to define the perturbation: abrupt or progressive, and which
	amplitude (+50%, +10.32%, or other).
        * lines 326-338: uncomment of one the lines to apply the perturbation in the desired case: additional carbonate dissolution,
	additional silicate dissolution, H2SO4 leaching.
	For the perturbation "silicate trade-off", ucomment the line 337 that deduces the additional flux from carbonic silicate
	weathering.
        * lines 427-434: for the perturbation "carbonate trade-off" only, uncomment the line 433 that deduces the additional flux
        from carbonic carbonate weathering.

    * Petrogenic carbon (kerogen) weathering perturbation:

        Still in the source file "cont_weath.f":
        * line 272-277: uncomment the line `fker(j) = ... * fker(j)` to apply the desired perturbation
        * line 391-395: if you want not to propagate the kerogen weathering perturbation to the phosphorus weathering flux,
	comment line 390 `+  P2C_ker * fker(j)`, and replace it by uncommenting one of the line `+  P2C_ker * fker(j) / ...`,
	where '...' is the kerogen weathering factor you applied.

	Both sulfide weathering and petrogenic carbon weathering perturbations can be applied simultaneously.

    * Oxygen feedback strength:
        * case "feedback-1": to impose a constant degree of anoxicity, comment the lines 6-7 in "DOfA.f" and uncomment line 15.
        * case "feedback-2": follow the instruction for "feedback-1", then in file "org_dep.f", comment line 55, and replace it
	by uncommenting line 57 ([O2]^0.5 dependency).
        * case no-feedback"": follow the instruction for "feedback-1", then in file "org_dep.f", comment line 55, and replace it
	by uncommenting line 59 (no O2 dependency).
        * case "feedback+1": add an [O2] dependency to hydrothermal P burial by commenting lines 8-9 in "Phydrotherm.f" and replacing
	them by uncommenting lines 11-12.
        * case "feedback+2": add an [O2]^2 dependency to hydrothermal P burial by commenting lines 8-9 in "Phydrotherm.f" and replacing
	them by uncommenting lines 14-15.
        * case "feedback+3": follow the instruction for "feedback+2", then in "cont_weath.f", comment line 290 and replace it by
	uncommenting line 292 (add a pO2^-0.5 dependency).
        * case "feedback+3": follow the instruction for "feedback+2", then in "cont_weath.f", comment line 290 and replace it by
	uncommenting line 294 (add a pO2^-1 dependency).

4. **Compile the code**:

    Once you have done the desired modifications in the source code, you can compile the code with `make` (since the pre-compilation
    configuration is already done with `configure.sh`), in the "source/" repertory:

    `make FC=your_fortran_compiler ncpath=path_to_netCDF_library`

    Alternatively, go back to the main repertory (`cd ..`) and compile with `build_GEOCLIM`:

    `./build_GEOCLIM --compset default --res 3,720,360`

    If you compile the code with `make`, the executable "geoclim.exe" will be created in the "source/" directory. I recommand to
    move it into the "executable/" directory.
    In any case, you should rename the executable just created to identify the code modifications you have made, especially if
    you want to do several modifications, and have to create as many executable as modifications. Otherwise, the executable will
    be overwritten at each new compilation.

5. **Chain run submission**:

    The most convenient way to run the simulations is to use the scripts for "chain" job submission, especially on a cluster, and
    if you wants to run several simulations in the same time.
    Those scripts ("submit_chain.sh" and "run_geoclim.sh") can be found in the "job/" directory

    * Go to the "job/" directory (`cd job`).
    * Replace the template "submit_chain.sh" by "submit_chain_sulfur-experiments.sh": `cp -f submit_chain_sulfur-experiments.sh
    submit_chain.sh`.  
    This last template contains the configuration for all the experiment presented in the article, whereas the first one is a
    default template.
    * Edit the "new" script "submit_chain.sh" to configure your run. The script is currently configure for the "carbonate sulfuric
    weathering" perturbation experiment. The configuration for all experiment are between the lines 38 and 238. You can simply
    comment the default config and uncomment the one corresponding to your run.
    The configuration variables that need to be edited are:
        * `RUN_NAME`: the name of your run
        * `EXECUTABLE`: the path & name or your GEOCLIM executable file 
        * `STOP_TIMES`, `COMBINE_DT`, `CONTWTH_NSKIP`, `DYNSOIL_NSKIP`, `COMBINE_PRINT_NSKIP`, `GEOGRAP_PRINT_NSKIP` and
	`DYNSOIL_PRINT_NSKIP`: variables (lists) controlling the times to stop and restart the runs, the timesteps of the different
	modules (COMBINE, continental weathering and DynSoil) and the printing timesteps.
	The default (uncommented) values are the ones that were used for all the abrupt perturbations.
        * `COMBINE_INIT`: the path & name for the initial condition of COMBINE variables. "restart/geoclim/output.ref" can be used
	for all simulations since they all start from a quasi identical pre-industrial steady-state.
	Yet, in order to start more closely to the numerical steady-state, different COMBINE restarts are available for all the
	oxygen feedback strength cases (see commented examples of configuration).
        * `DYNSOIL_INIT`: Path & name of the DynSoil initial condition file. Similarly to COMBINE initial condition, it would be
	preferrable to start from the numerical steady-state of the calibration run.
	However, DynSoil initialization files require too much memory to be stored on GitHub. For this reason, the current
	configuration starts from the analytical steady-state (computed by the code): `DYNSOIL_INIT='startup:eq'`.
	This inaccuracy of numerical versus analytical steady-state of the regolith profiles will cause a rapid (50kyr) peak of CO2
	of +15ppmv, resorbing in ~1Myr, if one launches a unperturbed run with 'startup:eq' DynSoil initial condition.  
	This "artificial" peak of CO2 is negligible with respect to the amplitude of the perturbations applied.
	To generate an exact steady-state DynSoil initial condition, re-compile and run the model without perturbation for 1Myr
        *with accelerated parameters* (i.e., S and O cycle accelerated x100, in "config/cond_p20.dat", and reduced regolith inertia:
	`scaling_factor` set to 1d-3 in "source/dynsoil\_physical\_parameters.f90")
        * `SUBMIT_COMMAND`: the command used to submit job on the cluster, if you are running the model on a cluster (on Cheyenne
	cluster, `qsub` (PBS) or `sbatch` (Slurm)).
	If you wants to directly run the model, without submitting a job, set `SUBMIT_COMMAND=''`.
        * `JOB_FILE`: the name of the 2nd script that actually run the executable. Let it 'run\_geoclim.sh' if you want to submit
	a job (you will likely need to edit the file). If you are not submitting a job (SUBMIT\_COMMAND=''), then, set it to
	'run\_geoclim\_basic.sh'
        * `LOG_FILE`: The name of the log file of the run. This is optional.
        * `GEOCLIM_IO_FILE`: The name of GEOCLIM main config file. Since all the configuration can be done "submit\_chain.sh",
	there is no need to use another file, so let it 'config/IO\_CONDITIONS'.
    * If you are running the model by submitting a job on a cluster, edit the script "run\_geoclim.sh". It should contain the "job"
    information (cluster account, required walltime, number of nodes wanted...). A template for Cheyenne cluster (using PBS) is
    available: "run\_geoclim\_cheyenne\_cluster.sh".  
    The script should normally contain the line `./geoclim.exe`, or `./geoclim.exe 0 1 3 0 0`. **Do not change the name**
    'geoclim.exe', it is an automatically-generated link toward the actual executable.  
    In any case, *the last line of the script must be* `test $? -eq 0 && ./submit_chain.sh`, that is the resubmission command.
    * Launch the run with `./submit_chain.sh`.  
    The script will automatically re-launch the successive runs. Several independent runs can be launched in parallel, as the
    script check for any usage conflict of the configuration files.
    If for some reason the submission is stalled (you keep receiving the message "configuration files busy. Run set-up postponed
    for 1 minute" while there is no other run launched), you can erase the security log with `./submit_chain.sh ERASE_CONFIG_QUEUE`.  
    The script "submit_chain.sh" creates a repertory named after the run name (usually, starting with '.', so it is a hidden
    repertory) that contains the main information of the chain of runs, *including the log file* (in the subrepertory "run/").
    Once one chain of runs is launched, you can safely modify the file "submit_chain.sh" to launch another serie of runs.
    * The outputs and restarts will be written in "OUTPUT/". The "intermediate" restarts will be automatically moved in
    "restart/geoclim/" and "restart/dynsoil/".
