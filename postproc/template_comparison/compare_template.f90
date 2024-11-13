program compare_template
!
! Load given run and compare it (for several key variables) to the corresponding template run.
! The run name is read in GEOCLIM's main file: "config/IO_CONDITIONS"

use netcdf
use netcdf_io_module, only: open_file, close_file, inquire_dim, nf90_check, get_var
use local_functions, only: compare_var1D, compare_var2D
implicit none

include 'path.inc' ! => get main path: 'geoclim_path'
include 'shape.inc' ! => get 'nbasin'
include 'output_size.inc' ! => get 'nCOMBoutvar'

! namelist variables
character(len=100):: run_name, check_run_name, phys_cond_file, killing_file
character(len=500):: output_directory, file_name
integer, parameter:: ndim=2
character(len=30):: vartype(nCOMBoutvar)
character(len=100):: dname(ndim), vname(nCOMBoutvar), units(nCOMBoutvar)
character(len=500):: long_name(nCOMBoutvar)
integer, dimension(nCOMBoutvar,ndim):: defdim
logical, dimension(nCOMBoutvar):: writevar
double precision, dimension(nCOMBoutvar):: fillval

! main variables
character(len=1000):: fname, fname_0
integer, parameter:: nvar=24
real, dimension(:), allocatable:: CO2_level,   O2_level,   sil_wth,   bas_wth,   carb_wth,   ker_wth,   phos_wth, &
                                  terr_bio_C_exp,   pel_carb_dep,   time
real, dimension(:), allocatable:: CO2_level_0, O2_level_0, sil_wth_0, bas_wth_0, carb_wth_0, ker_wth_0, phos_wth_0, &
                                  terr_bio_C_exp_0, pel_carb_dep_0, time_0
real, dimension(:,:), allocatable:: O2,   DIC,   alk,   pH,   Sr,   Sr_ratio,   d13C,   bioprod,   ner_carb_dep,   sedim_flux, &
                                    org_C_dep,   bur_eff,   P_org_dep,   P_phos_dep,   P_hydro_dep
real, dimension(:,:), allocatable:: O2_0, DIC_0, alk_0, pH_0, Sr_0, Sr_ratio_0, d13C_0, bioprod_0, ner_carb_dep_0, sedim_flux_0, &
                                    org_C_dep_0, bur_eff_0, P_org_dep_0, P_phos_dep_0, P_hydro_dep_0
real:: delta, xref
integer:: fid, fid_0, dimid, ierr
integer:: nt, nt0, k, j
logical, dimension(nvar+1):: passed

! Namelist declaration
namelist /MAIN_INFO/ run_name, output_directory, phys_cond_file, killing_file
namelist /CMB_OUTPUT_FILE/ file_name
namelist /CMB_OUTPUT_DIM/ dname, units, vartype
namelist /CMB_OUTPUT_VAR/ vname, units, defdim, writevar, long_name, fillval, vartype




print *
print *, 'GEOCLIM template-reference comparison check'
print *, '###########################################'
print *
print *




! Load GEOCLIM OUTPUT information
! ===============================


! Open IO_CONDITION file
open(unit=1, status='old', action='read', file=geoclim_path//'config/IO_CONDITIONS')

! Read namelists
read(unit=1, nml=MAIN_INFO)
read(unit=1, nml=CMB_OUTPUT_FILE)
read(unit=1, nml=CMB_OUTPUT_DIM)
read(unit=1, nml=CMB_OUTPUT_VAR)

! close IO_CONDITION file
close(unit=1)



! Load variables
! ==============


! TIME
!
!     - new output file
fname = geoclim_path//'OUTPUT/'//trim(file_name)//trim(run_name)//'.nc'
call open_file(fname, fid)
call inquire_dim(fid, trim(dname(2)), dimid)
ierr = nf90_inquire_dimension(fid, dimid, len=nt)
call nf90_check(ierr, 'Error while inquiring ID of variable "'//trim(dname(2))//'" in file "'//trim(fname)//'"')
!
!     - Reference file:
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(file_name)//trim(run_name)//'.nc'
call open_file(fname_0, fid_0)
call inquire_dim(fid_0, trim(dname(2)), dimid)
ierr = nf90_inquire_dimension(fid_0, dimid, len=nt0)
call nf90_check(ierr, 'Error while inquiring ID of variable "'//trim(dname(2))//'" in file "'//trim(fname_0)//'"')
!
! Compare length
if (nt/=nt0) then
    print *
    print *, 'TEST NOT PASSED: time dimensions do not have the same length'
    stop
end if
!
!################################!
! ALLOCATE VARIABLES
allocate(time(nt))
allocate(time_0(nt))
allocate(CO2_level(nt))
allocate(CO2_level_0(nt))
allocate(O2_level(nt))
allocate(O2_level_0(nt))
allocate(sil_wth(nt))
allocate(sil_wth_0(nt))
allocate(bas_wth(nt))
allocate(bas_wth_0(nt))
allocate(carb_wth(nt))
allocate(carb_wth_0(nt))
allocate(ker_wth(nt))
allocate(ker_wth_0(nt))
allocate(phos_wth(nt))
allocate(phos_wth_0(nt))
allocate(terr_bio_C_exp(nt))
allocate(terr_bio_C_exp_0(nt))
allocate(pel_carb_dep(nt))
allocate(pel_carb_dep_0(nt))
allocate(O2(nbasin,nt))
allocate(O2_0(nbasin,nt))
allocate(DIC(nbasin,nt))
allocate(DIC_0(nbasin,nt))
allocate(alk(nbasin,nt))
allocate(alk_0(nbasin,nt))
allocate(pH(nbasin,nt))
allocate(pH_0(nbasin,nt))
allocate(Sr(nbasin,nt))
allocate(Sr_0(nbasin,nt))
allocate(Sr_ratio(nbasin,nt))
allocate(Sr_ratio_0(nbasin,nt))
allocate(d13C(nbasin,nt))
allocate(d13C_0(nbasin,nt))
allocate(bioprod(nbasin,nt))
allocate(bioprod_0(nbasin,nt))
allocate(ner_carb_dep(nbasin,nt))
allocate(ner_carb_dep_0(nbasin,nt))
allocate(sedim_flux(nbasin,nt))
allocate(sedim_flux_0(nbasin,nt))
allocate(org_C_dep(nbasin,nt))
allocate(org_C_dep_0(nbasin,nt))
allocate(bur_eff(nbasin,nt))
allocate(bur_eff_0(nbasin,nt))
allocate(P_org_dep(nbasin,nt))
allocate(P_org_dep_0(nbasin,nt))
allocate(P_phos_dep(nbasin,nt))
allocate(P_phos_dep_0(nbasin,nt))
allocate(P_hydro_dep(nbasin,nt))
allocate(P_hydro_dep_0(nbasin,nt))
!################################!
!
! load time variable
call get_var(fid,   dname(2), var_real1D=time)
call get_var(fid_0, dname(2), var_real1D=time_0)


! ATMOSPHERIC CO2
call get_var(fid,   vname(58), var_real1D=CO2_level)
call get_var(fid_0, vname(58), var_real1D=CO2_level_0)

! ATMOSPHERIC O2
call get_var(fid,   vname(57), var_real1D=O2_level)
call get_var(fid_0, vname(57), var_real1D=O2_level_0)

! SILICATE WEATHERING
call get_var(fid,   vname(67), var_real1D=sil_wth)
call get_var(fid_0, vname(67), var_real1D=sil_wth_0)

! BASALT WEATHERINg
call get_var(fid,   vname(68), var_real1D=bas_wth)
call get_var(fid_0, vname(68), var_real1D=bas_wth_0)

! CARBONATE WEATHERING
call get_var(fid,   vname(69), var_real1D=carb_wth)
call get_var(fid_0, vname(69), var_real1D=carb_wth_0)

! KEROGEN WEATHERING
call get_var(fid,   vname(70), var_real1D=ker_wth)
call get_var(fid_0, vname(70), var_real1D=ker_wth_0)

! PHOSPHORUS WEATHERING
call get_var(fid,   vname(77), var_real1D=phos_wth)
call get_var(fid_0, vname(77), var_real1D=phos_wth_0)

! TERRESTRIAL BIOGENIC C EXPORTS
call get_var(fid,   vname(87), var_real1D=terr_bio_C_exp)
call get_var(fid_0, vname(87), var_real1D=terr_bio_C_exp_0)

! PELAGIC CARBONATE DEPOSITION
call get_var(fid,   vname(73), var_real1D=pel_carb_dep)
call get_var(fid_0, vname(73), var_real1D=pel_carb_dep_0)

! DISSOLVED O2
call get_var(fid,   vname(14), var_real2D=O2)
call get_var(fid_0, vname(14), var_real2D=O2_0)

! DIC
call get_var(fid,   vname(4), var_real2D=DIC)
call get_var(fid_0, vname(4), var_real2D=DIC_0)

! ALKALINITY
call get_var(fid,   vname(5), var_real2D=alk)
call get_var(fid_0, vname(5), var_real2D=alk_0)

! PH
call get_var(fid,   vname(30), var_real2D=pH)
call get_var(fid_0, vname(30), var_real2D=pH_0)

! SR CONCENTRATION
call get_var(fid,   vname(8), var_real2D=Sr)
call get_var(fid_0, vname(8), var_real2D=Sr_0)

! SR ISOTOPIC RATIO
call get_var(fid,   vname(20), var_real2D=Sr_ratio)
call get_var(fid_0, vname(20), var_real2D=Sr_ratio_0)

! DELTA 13 C
call get_var(fid,   vname(16), var_real2D=d13C)
call get_var(fid_0, vname(16), var_real2D=d13C_0)

! BIOPRODUCTIVITY
call get_var(fid,   vname(78), var_real2D=bioprod)
call get_var(fid_0, vname(78), var_real2D=bioprod_0)

! NERITIC CARBONATE DEPOPSITION
call get_var(fid,   vname(71), var_real2D=ner_carb_dep)
call get_var(fid_0, vname(71), var_real2D=ner_carb_dep_0)

! SEDIMENTATION FLUX
call get_var(fid,   vname(95), var_real2D=sedim_flux)
call get_var(fid_0, vname(95), var_real2D=sedim_flux_0)

! ORGANIC C DEPOSITION
call get_var(fid,   vname(74), var_real2D=org_C_dep)
call get_var(fid_0, vname(74), var_real2D=org_C_dep_0)

! BURIAL EFFICIENCY
call get_var(fid,   vname(96), var_real2D=bur_eff)
call get_var(fid_0, vname(96), var_real2D=bur_eff_0)

! ORGANIC P DEPOSITION
call get_var(fid,   vname(97), var_real2D=P_org_dep)
call get_var(fid_0, vname(97), var_real2D=P_org_dep_0)

! PHOSPHORITE P DEPOSITION
call get_var(fid,   vname(98), var_real2D=P_phos_dep)
call get_var(fid_0, vname(98), var_real2D=P_phos_dep_0)

! HYDROTHERMAL P DEPOSITION
call get_var(fid,   vname(99), var_real2D=P_hydro_dep)
call get_var(fid_0, vname(99), var_real2D=P_hydro_dep_0)


! Close netCDF files
call close_file(fid)
call close_file(fid_0)




! Compare variables and print report
! ==================================

print *
print *


passed(1)  = compare_var1D('time',                 time,           time_0)
passed(2)  = compare_var1D('atmospheric CO2',      CO2_level,      CO2_level_0)
passed(3)  = compare_var1D('atmospheric O2',       O2_level,       O2_level_0)
passed(4)  = compare_var1D('silicate weathering',  sil_wth,        sil_wth_0)
passed(5)  = compare_var1D('basalt weathering',    bas_wth,        bas_wth_0)
passed(6)  = compare_var1D('carbonate weathering', carb_wth,       carb_wth_0)
passed(7)  = compare_var1D('kerogen weathering',   ker_wth,        ker_wth_0)
passed(8)  = compare_var1D('phosphorus weather.',  phos_wth,       phos_wth_0)
passed(9)  = compare_var1D('terr. bio. C export',  terr_bio_C_exp, terr_bio_C_exp_0)
passed(10) = compare_var1D('pelagic carb. depos.', pel_carb_dep,   pel_carb_dep_0)
passed(11) = compare_var2D('dissolved O2',         O2,             O2_0)
passed(12) = compare_var2D('DIC',                  DIC,            DIC_0)
passed(13) = compare_var2D('Alkalinity',           alk,            alk_0)
passed(14) = compare_var2D('pH',                   pH,             pH_0)
passed(15) = compare_var2D('Sr concentration',     Sr,             Sr_0)
passed(16) = compare_var2D('Sr isotopic ratio',    Sr_ratio,       Sr_ratio_0)
passed(17) = compare_var2D('DIC delta 13 C',       d13C,           d13C_0)
passed(18) = compare_var2D('bioproductivity',      bioprod,        bioprod_0)
passed(19) = compare_var2D('neritic carb. depos.', ner_carb_dep,   ner_carb_dep_0)
passed(20) = compare_var2D('sedimentation flux',   sedim_flux,     sedim_flux_0)
passed(21) = compare_var2D('organic C deposition', org_C_dep,      org_C_dep_0)
passed(22) = compare_var2D('C burial efficiency',  bur_eff,        bur_eff_0)
passed(23) = compare_var2D('organic P deposition', P_org_dep,      P_org_dep_0)
passed(24) = compare_var2D('phosphorite P depos.', P_phos_dep,     P_phos_dep_0)
passed(25) = compare_var2D('hydrothermal P dep.',  P_hydro_dep,    P_hydro_dep_0)




print *
print *
if (all(passed)) then
    print *, 'Template comparison check PASSED'
else
    print *, 'Template comparison check NOT PASSED'
end if



end program
