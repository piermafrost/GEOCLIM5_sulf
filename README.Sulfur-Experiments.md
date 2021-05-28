# Sulfide weathering experiments from Maffre et al., submitted to GRL (2021)

Information about GEOCLIM runs for the sulfide weathering perturbation experiments

## Run outputs

#### Where to find the output

The outputs of the sulfide weathering simulations are stored in "OUTPUT/sulfur\_experiments/".
Only the COMBINE outputs (i.e., ocean and atmosphere "box" variables and integrated fluxes) are available.
The geographic outputs are too heavy be be stored in GitHub.

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

    `./build_GEOCLIM --compset default --res 3,360,720`

    If you compile the code with `make`, the executable "geoclim.exe" will be created in the "source/" directory. I recommand to
    move it into the "executable/" directory.
    In any case, you should rename the executable just created to identify the code modifications you have made, especially if
    you want to do several modifications, and have to create as many executable as modifications. Otherwise, the executable will
    be overwritten at each new compilation.

5. **chain job submission**:

    Hey?
