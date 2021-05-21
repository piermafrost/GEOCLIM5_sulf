module local_functions
implicit none

! Maximum allowed relative error
real, parameter:: MAX_ALLOWED_ERROR = 1e-6


  contains



    function compare_var1D(varname, var, var_ref)

        logical:: compare_var1D
        character(len=*), intent(in):: varname
        real, dimension(:), intent(in):: var, var_ref
        logical:: passed
        real:: xref, delta
        integer:: nt

        nt = size(var)

        xref = abs(sum(var_ref)) / nt
        delta = maxval(abs(var - var_ref))
        if (xref > 1e-12) delta = delta/xref
        compare_var1D = (delta <= MAX_ALLOWED_ERROR)

        ! Print message
        if (compare_var1D) then
            write(*, fmt='(A4,A40,ES10.3,A1)') '  * ', trim(varname)//': PASSED   (delta = ', delta, ')'
        else
            write(*, fmt='(A4,A40,ES10.3)')     '  * ', trim(varname)//': FAILED:   delta = ', delta
        end if
        write(*, fmt='(A)') '                          ----------'

    end function


    !--------!


    function compare_var2D(varname, var, var_ref, sedim)

        include 'shape.inc' ! => 'nbasin'
        integer, dimension(5), parameter:: sedim_basins = (/2, 5, 6, 7, 9/)

        logical:: compare_var2D
        character(len=*), intent(in):: varname
        real, dimension(:,:), intent(in):: var, var_ref
        logical, intent(in), optional:: sedim
        logical:: loc_sedim, passed
        real:: xref, loc_delta
        real, dimension(nbasin):: delta
        integer, dimension(nbasin):: failed
        integer:: j, j0, jend, nt, k

        nt = size(var, 2)
        k = 0

        if (present(sedim)) then
            loc_sedim = sedim
        else
            loc_sedim = .false.
        end if

        if (loc_sedim) then
            jend = size(sedim_basins)
            do j0 = 1,jend
                j = sedim_basins(j0)
                xref = abs(sum(var_ref(j,:))) / nt
                loc_delta = maxval(abs(var(j,:) - var_ref(j,:)))
                if (xref > 1e-12) loc_delta = loc_delta/xref
                passed = (loc_delta <= MAX_ALLOWED_ERROR)
                if (.not. passed) then
                    k = k + 1
                    failed(k) = j
                    delta(k) = loc_delta
                end if
            end do
        else
            jend = nbasin-1 ! skip atmospheric box
            do j = 1,jend
                xref = abs(sum(var_ref(j,:))) / nt
                loc_delta = maxval(abs(var(j,:) - var_ref(j,:)))
                if (xref > 1e-12) loc_delta = loc_delta/xref
                passed = (loc_delta <= MAX_ALLOWED_ERROR)
                if (.not. passed) then
                    k = k + 1
                    failed(k) = j
                    delta(k) = loc_delta
                end if
            end do
        end if

        ! Print message
        if (k==0) then
            write(*, fmt='(A4,A40,ES10.3,A1)') '  * ', trim(varname)//': PASSED   (delta = ', maxval(delta(1:jend)), ')'
            compare_var2D = .true.
        else
            write(*, fmt='(A4,A28,ES10.3)') '  * ', trim(varname)//': FAILED'
            do j = 1,k
                write(*, fmt='(A14,I0,A9,ES10.3)') '       basin #', failed(j), 'delta = ', delta(j)
            end do
            compare_var2D = .false.
        end if
        write(*, fmt='(A)') '                          ----------'

    end function


end module





program compare_template
!
! Load given run and compare it (for several key variables) to the corresponding template run.
! The run name is read in GEOCLIM's main file: "config/IO_CONDITIONS"

use netcdf
use utils, only: read_comment
use netcdf_io_module, only: open_file, close_file, inquire_dim, isitok, get_var_real1D, get_var_real2D
use local_functions, only: compare_var1D, compare_var2D
implicit none

include 'path.inc' ! => get main path: 'geoclim_path'
include 'shape.inc' ! => get 'nbasin'


character(len=1000):: fname, fname_0
character(len=200):: run_name, time_varname, dummy
integer, parameter:: nvar=24
character(len=200), dimension(nvar):: output_file, varname
real, dimension(:), allocatable:: CO2_level,   O2_level,   sil_wth,   bas_wth,   carb_wth,   ker_wth,   phos_wth, &
                                  terr_bio_C_exp,   pel_carb_dep,   time
real, dimension(:), allocatable:: CO2_level_0, O2_level_0, sil_wth_0, bas_wth_0, carb_wth_0, ker_wth_0, phos_wth_0, &
                                  terr_bio_C_exp_0, pel_carb_dep_0, time_0
real, dimension(:,:), allocatable:: O2,   DIC,   alk,   pH,   Sr,   Sr_ratio,   d13C,   bioprod,   ner_carb_dep,   sedim_flux, &
                                    org_C_dep,   bur_eff,   P_org_dep,   P_phos_dep,   P_hydro_dep
real, dimension(:,:), allocatable:: O2_0, DIC_0, alk_0, pH_0, Sr_0, Sr_ratio_0, d13C_0, bioprod_0, ner_carb_dep_0, sedim_flux_0, &
                                    org_C_dep_0, bur_eff_0, P_org_dep_0, P_phos_dep_0, P_hydro_dep_0
real:: delta, xref
integer:: fid, fid_0, dimid(1), ierr
integer:: nt, nt0, k, j
logical, dimension(nvar+1):: passed




print *
print *, 'GEOCLIM template-reference comparison check'
print *, '###########################################'
print *
print *




! Load GEOCLIM OUTPUT information
! ===============================


! Open IO_CONDITION file
open(unit=1, status='old', action='read', file=geoclim_path//'config/IO_CONDITIONS')

! Read run name
call read_comment(1)
read(unit=1, fmt=*) dummy, run_name

! Read variables name and output file:

! * time
do k=1,58; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, dummy, time_varname
!
! * DIC
do k=1,3; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(11), varname(11)
!
! * Alkalinity
read(unit=1, fmt=*) dummy, output_file(12), varname(12)
!
! * Strontium
do k=1,2; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(14), varname(14)
!
! * O2
do k=1,5; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(10), varname(10)
!
! * d13C
call read_comment(1); read(unit=1, fmt=*)
read(unit=1, fmt=*) dummy, output_file(16), varname(16)
!
! * Strontium ratio
do k=1,3; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(15), varname(15)
!
! * pH
do k=1,9; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(13), varname(13)
!
! * O2 atmospheric level
do k=1,27; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(2), varname(2)
!
! * CO2 atmospheric level
read(unit=1, fmt=*) dummy, output_file(1), varname(1)
!
! * Silicate weathering
do k=1,8; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(3), varname(3)
!
! * Basalt weathering
read(unit=1, fmt=*) dummy, output_file(4), varname(4)
!
! * Carbonate weathering
read(unit=1, fmt=*) dummy, output_file(5), varname(5)
!
! * Kerogen weathering
read(unit=1, fmt=*) dummy, output_file(6), varname(6)
!
! * Neritic carbonae deposition
read(unit=1, fmt=*) dummy, output_file(18), varname(18)
!
! * Total pelagic carbonate deposition 
call read_comment(1); read(unit=1, fmt=*)
read(unit=1, fmt=*) dummy, output_file(9), varname(9)
!
! * Organic carbon burial
read(unit=1, fmt=*) dummy, output_file(20), varname(20)
!
! * Phosphorus weathering
do k=1,2; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(7), varname(7)
!
! * Bioproductivity
read(unit=1, fmt=*) dummy, output_file(17), varname(17)
!
! * Terrestrial biogenic POC export
do k=1,8; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(8), varname(8)
!
! * Sedimentation flux
do k=1,7; call read_comment(1); read(unit=1, fmt=*); end do
read(unit=1, fmt=*) dummy, output_file(19), varname(19)
!
! * Organic carbon burial efficiency
read(unit=1, fmt=*) dummy, output_file(21), varname(21)
!
! * Organic carbon-bound phosphorus burial
read(unit=1, fmt=*) dummy, output_file(22), varname(22)
!
! * Phosphorus burial in form of phophorite
read(unit=1, fmt=*) dummy, output_file(23), varname(23)
!
! * Hydrothermal Fe-bound Phosphorus burial
read(unit=1, fmt=*) dummy, output_file(24), varname(24)


! close IO_CONDITION file
close(unit=1)


! 1st file (main) is DIC (loc #11):
if (output_file(1)=='-') output_file(1) = output_file(11)

do k=2,24
    ! replace "-" by previous line's file name
    dummy = output_file(k)
    j = k
    do while (dummy=='-')
        j = j-1
        dummy = output_file(j)
    end do
    if (dummy == '#') then
        print *
        print *, 'ERROR: missing file for output variable '//trim(varname(j))
        print *, 'Cannot make comparison'
        stop
    else
        output_file(k) = output_file(j)
    end if
end do




! Load variables
! ==============


! TIME
!
!     - new output file
fname = geoclim_path//'OUTPUT/'//trim(output_file(1))//trim(run_name)//'.nc'
call open_file(fname, fid)
call inquire_dim(fid, (/trim(time_varname)/), dimid)
ierr = nf90_inquire_dimension(fid, dimid(1), len=nt)
call isitok(ierr, 'Error while inquiring ID of variable "'//trim(time_varname)//'" in file "'//trim(fname)//'"')
!
!     - Reference file:
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(1))//trim(run_name)//'.nc'
call open_file(fname_0, fid_0)
call inquire_dim(fid_0, (/trim(time_varname)/), dimid)
ierr = nf90_inquire_dimension(fid_0, dimid(1), len=nt0)
call isitok(ierr, 'Error while inquiring ID of variable "'//trim(time_varname)//'" in file "'//trim(fname_0)//'"')
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
! load variable
call get_var_real1D(fid, time_varname, time)
call get_var_real1D(fid_0, time_varname, time_0)
call close_file(fid)
call close_file(fid_0)


! ATMOSPHERIC CO2
fname   = geoclim_path//'OUTPUT/'//trim(output_file(1))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(1))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(1), CO2_level)
call get_var_real1D(fid_0, varname(1), CO2_level_0)
call close_file(fid)
call close_file(fid_0)

! ATMOSPHERIC O2
fname   = geoclim_path//'OUTPUT/'//trim(output_file(2))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(2))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(2), O2_level)
call get_var_real1D(fid_0, varname(2), O2_level_0)
call close_file(fid)
call close_file(fid_0)

! SILICATE WEATHERING
fname   = geoclim_path//'OUTPUT/'//trim(output_file(3))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(3))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(3), sil_wth)
call get_var_real1D(fid_0, varname(3), sil_wth_0)
call close_file(fid)
call close_file(fid_0)

! BASALT WEATHERINg
fname   = geoclim_path//'OUTPUT/'//trim(output_file(4))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(4))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(4), bas_wth)
call get_var_real1D(fid_0, varname(4), bas_wth_0)
call close_file(fid)
call close_file(fid_0)

! CARBONATE WEATHERING
fname   = geoclim_path//'OUTPUT/'//trim(output_file(5))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(5))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(5), carb_wth)
call get_var_real1D(fid_0, varname(5), carb_wth_0)
call close_file(fid)
call close_file(fid_0)

! KEROGEN WEATHERING
fname   = geoclim_path//'OUTPUT/'//trim(output_file(6))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(6))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(6), ker_wth)
call get_var_real1D(fid_0, varname(6), ker_wth_0)
call close_file(fid)
call close_file(fid_0)

! PHOSPHORUS WEATHERIN
fname   = geoclim_path//'OUTPUT/'//trim(output_file(7))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(7))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(7), phos_wth)
call get_var_real1D(fid_0, varname(7), phos_wth_0)
call close_file(fid)
call close_file(fid_0)

! TERRESTRIAL BIOGENIC C EXPORTS
fname   = geoclim_path//'OUTPUT/'//trim(output_file(8))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(8))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(8), terr_bio_C_exp)
call get_var_real1D(fid_0, varname(8), terr_bio_C_exp_0)
call close_file(fid)
call close_file(fid_0)

! PELAGIC CARBONATE DEPOSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(9))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(9))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real1D(fid,   varname(9), pel_carb_dep)
call get_var_real1D(fid_0, varname(9), pel_carb_dep_0)
call close_file(fid)
call close_file(fid_0)

! DISSOLVED O2
fname   = geoclim_path//'OUTPUT/'//trim(output_file(10))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(10))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(10), O2)
call get_var_real2D(fid_0, varname(10), O2_0)
call close_file(fid)
call close_file(fid_0)

! DIC
fname   = geoclim_path//'OUTPUT/'//trim(output_file(11))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(11))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(11), DIC)
call get_var_real2D(fid_0, varname(11), DIC_0)
call close_file(fid)
call close_file(fid_0)

! ALKALINITY
fname   = geoclim_path//'OUTPUT/'//trim(output_file(12))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(12))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(12), alk)
call get_var_real2D(fid_0, varname(12), alk_0)
call close_file(fid)
call close_file(fid_0)

! PH
fname   = geoclim_path//'OUTPUT/'//trim(output_file(13))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(13))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(13), pH)
call get_var_real2D(fid_0, varname(13), pH_0)
call close_file(fid)
call close_file(fid_0)

! SR CONCENTRATION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(14))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(14))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(14), Sr)
call get_var_real2D(fid_0, varname(14), Sr_0)
call close_file(fid)
call close_file(fid_0)

! SR ISOTOPIC RATIO
fname   = geoclim_path//'OUTPUT/'//trim(output_file(15))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(15))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(15), Sr_ratio)
call get_var_real2D(fid_0, varname(15), Sr_ratio_0)
call close_file(fid)
call close_file(fid_0)

! DELTA 13 C
fname   = geoclim_path//'OUTPUT/'//trim(output_file(16))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(16))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(16), d13C)
call get_var_real2D(fid_0, varname(16), d13C_0)
call close_file(fid)
call close_file(fid_0)

! BIOPRODUCTIVITY
fname   = geoclim_path//'OUTPUT/'//trim(output_file(17))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(17))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(17), bioprod)
call get_var_real2D(fid_0, varname(17), bioprod_0)
call close_file(fid)
call close_file(fid_0)

! NERITIC CARBONATE DEPOPSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(18))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(18))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(18), ner_carb_dep)
call get_var_real2D(fid_0, varname(18), ner_carb_dep_0)
call close_file(fid)
call close_file(fid_0)

! SEDIMENTATION FLUX
fname   = geoclim_path//'OUTPUT/'//trim(output_file(19))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(19))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(19), sedim_flux)
call get_var_real2D(fid_0, varname(19), sedim_flux_0)
call close_file(fid)
call close_file(fid_0)

! ORGANIC C DEPOSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(20))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(20))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(20), org_C_dep)
call get_var_real2D(fid_0, varname(20), org_C_dep_0)
call close_file(fid)
call close_file(fid_0)

! BURIAL EFFICIENCY
fname   = geoclim_path//'OUTPUT/'//trim(output_file(21))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(21))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(21), bur_eff)
call get_var_real2D(fid_0, varname(21), bur_eff_0)
call close_file(fid)
call close_file(fid_0)

! ORGANIC P DEPOSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(22))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(22))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(22), P_org_dep)
call get_var_real2D(fid_0, varname(22), P_org_dep_0)
call close_file(fid)
call close_file(fid_0)

! PHOSPHORITE P DEPOSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(23))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(23))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(23), P_phos_dep)
call get_var_real2D(fid_0, varname(23), P_phos_dep_0)
call close_file(fid)
call close_file(fid_0)

! HYDROTHERMAL P DEPOSITION
fname   = geoclim_path//'OUTPUT/'//trim(output_file(24))//trim(run_name)//'.nc'
fname_0 = geoclim_path//'OUTPUT/templates/'//trim(output_file(24))//trim(run_name)//'.nc'
call open_file(fname,   fid)
call open_file(fname_0, fid_0)
call get_var_real2D(fid,   varname(24), P_hydro_dep)
call get_var_real2D(fid_0, varname(24), P_hydro_dep_0)
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
passed(15) = .true. !compare_var2D('Sr concentration',     Sr,             Sr_0)
passed(16) = .true. !compare_var2D('Sr isotopic ratio',    Sr_ratio,       Sr_ratio_0)
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
