module geoclim_write_output_mod
implicit none

contains

subroutine geoclim_write_output(t, COMB_outvar_info)
  use netcdf
  use io_module,        only: netcdf_output_var
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, put_var
  use constante,        only: gotoshelf, rho_sed, PI_n_O2_atm, PI_n_CO2_atm, PI_n_rest_of_atm
  include 'combine_foam.inc'
  integer:: ierr, fid, timevarid, dimid, nt
  ! Declare structures for netCDF output info (cannot be declared in combine_foam.inc)
  type(netcdf_output_var), dimension(nCOMBoutvar), intent(in) :: COMB_outvar_info


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!           offline computations:           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xPOPexport=0.
    do j0=1,nsurface
        j = jbox_surface(j0)
        xPOPexport=xPOPexport + fsink(j)*var_part(2,j)*box_vol(j)
        !totalCO2_exch=totalCO2_exch+fco2atm_ocean(j)
    enddo

    vol_tot=0.
    vartot_diss(:)=0.
    vartot_part(:)=0.
    vartot_isot(:)=0.
    ph_tot=0.
    temp_tot=0.
    salin_tot=0.
    freef_tot=0.
    fodc_tot=0.
    fdep_tot=0.
    fco2atm_ocean_tot=0.
    F_seafloor_cdiss_tot=0.
    fSulfRed_tot=0.
    fO2_odc_tot=0.
    do j=1,nbasin-1
        vol_tot              = vol_tot + box_vol(j)
        vartot_diss(:)       = vartot_diss(:) + var_diss(:,j)*box_vol(j) ! dissolved variables
        vartot_part(:)       = vartot_part(:) + var_part(:,j)*box_vol(j) ! particulate variables
        vartot_isot(1)       = vartot_isot(1) + var_isot(1,j)*var_diss(1,j)*box_vol(j) ! DIC d13C
        vartot_isot(2)       = vartot_isot(2) + var_isot(2,j)*var_part(5,j)*box_vol(j) ! PIC d13C
        vartot_isot(3)       = vartot_isot(3) + var_isot(3,j)*var_part(4,j)*box_vol(j) ! POC d13C
        vartot_isot(5)       = vartot_isot(5) + var_isot(5,j)*var_diss(5,j)*box_vol(j) ! diss Sr isotopic ration
        vartot_isot(6)       = vartot_isot(6) + var_isot(6,j)*var_part(1,j)*box_vol(j) ! PIC Sr isotopic ration
        !vartot_isot()       = vartot_isot() + var_isot(,j)*var_diss(,j)*box_vol(j) ! d7Li
        ph_tot               = ph_tot + pH(j)*box_vol(j)
        temp_tot             = temp_tot + temp_box(j)*box_vol(j)
        salin_tot            = salin_tot + salin(j)*box_vol(j)
        freef_tot            = freef_tot+freef(j)
        fodc_tot             = fodc_tot+fodc(j)
        fco2atm_ocean_tot    = fco2atm_ocean_tot + fco2atm_ocean(j)
        F_seafloor_cdiss_tot = F_seafloor_cdiss_tot + F_seafloor_cdiss(j)
        fSulfRed_tot         = fSulfRed_tot + fSulfRed(j)
        fO2_odc_tot          = fO2_odc_tot + fO2_odc(j)
    end do
    ! Carbonate deposition: sum "fsink*[PIC]" on sedimentary reservoirs
    ! + remove epicont surf fraction of fsink that did not go into shelf (1-gotoshelf)
    do j0=1,nsedi
        j = jbox_sedi(j0)
        fdep_tot = fdep_tot + fsink_inorg(j)*var_part(5,j)*(box_vol(j))
    end do
    do j0=1,nappcont-1!<- skip atmospheric box, counted in "appcont"
        j = jbox_appcont(j0)
        fdep_tot = fdep_tot - (1-gotoshelf)*fsink_inorg(j)*var_part(5,j)*(box_vol(j))
    end do

    ! Normalize isotopic variables
    vartot_isot(1) = vartot_isot(1) / vartot_diss(1) ! DIC d13C
    !
    if (vartot_part(5) == 0d0) then
        vartot_isot(2) = COMB_outvar_info(14)%fillval
    else
        vartot_isot(2) = vartot_isot(2) / vartot_diss(5) ! PIC d13C
    end if
    !
    vartot_isot(3) = vartot_isot(3) / vartot_part(4) ! POC d13C
    !
    if (vartot_diss(5) == 0d0) then
        vartot_isot(5) = COMB_outvar_info(17)%fillval
    else
        vartot_isot(5) = vartot_isot(5) / vartot_diss(5) ! diss Sr iso ration
    end if
    !
    if (vartot_part(1) == 0d0) then
        vartot_isot(6) = COMB_outvar_info(18)%fillval
    else
        vartot_isot(6) = vartot_isot(6) / vartot_part(1) ! PIC Sr iso ration
    end if
    !
    !vartot_isot()    = vartot_isot() / vartot_diss() ! d7Li

    ! Normalize mean intensive variables (e.g., concentration)
    vartot_diss    = vartot_diss / vol_tot ! dissolved variables
    vartot_part    = vartot_part / vol_tot ! particulate variables

    ph_tot         = ph_tot / vol_tot
    temp_tot       = temp_tot / vol_tot
    salin_tot      = salin_tot / vol_tot

    ! Present-day normalized mol of O2 and CO2
    O2_atm_level  = var_diss(6,nbasin) / PI_n_O2_atm
    CO2_atm_level = var_diss(7,nbasin) / PI_n_CO2_atm

    ! convert mol to volume fraction (= molar fraction) for 02 and CO2
    O2_atm_conc   =  ( var_diss(6,nbasin) / ( PI_n_rest_of_atm + var_diss(6,nbasin) + var_diss(7,nbasin) ) )*100 ! (%)
    CO2_atm_conc  =  ( var_diss(7,nbasin) / ( PI_n_rest_of_atm + var_diss(6,nbasin) + var_diss(7,nbasin) ) )*1e6 ! (ppm)



!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!


    ! Open file, get size of time dimension, and put current time
    ! -----------------------------------------------------------
    call open_file(COMB_ofile_name, fid, mode=NF90_WRITE)
    call inquire_dim(fid, COMB_time_dimname, dimid)
    ierr = nf90_inquire_dimension(fid, dimid, len=nt)
    nt = nt + 1
    call inquire_var(fid, COMB_time_dimname, timevarid)
    call put_var(fid, timevarid, var_real0D=real(t), stt=(/nt/), cnt=(/1/))



    !<><><><><><><><><><><><><><><><><>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<> %  write output variables  % <>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<><><><><><><><><><><><><><><><><>!

    ! Note: the following list of blocks needs to be updated if one wants to add new output variables
    ! ***********************************************************************************************

    do i = 4,8 ! => first 5 main COMBINE dissolved variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_diss(i-3,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    do i = 9,13 ! => 5 main COMBINE particulate variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_part(i-8,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    do i = 14,15 ! => main COMBINE dissolved variables # 6 and 7 (O2 and pCO2)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_diss(i-8,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    do i = 16,21 ! => 6 main COMBINE isotopic variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_isot(i-15,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
    !
    ! i=22 corresponds to former unattributed variable #19
    !
    i = 23 ! main COMBINE variable # 8 (sulfate)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var_diss(8,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 24
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(h2co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 25
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(hco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 26
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 27
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dh2co3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 28
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dhco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 29
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dco3(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 30
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(ph(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 31
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(omega(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 32
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(temp_box(1:nbasin-1)-273.15), &
                     stt=(/1,nt/), cnt=(/nbasin-1,1/))
    !
    i = 33
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(salin(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 34
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dplysc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 35
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(dplysa(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    do i = 36,40 ! => box-average values of the main COMBINE dissolved variables #1 to #5
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_diss(i-35)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 41,45 ! => box-average values of the 5 main COMBINE particulate variables
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_part(i-40)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 46 ! => box-average value of the main COMBINE dissolved variable #6 (O2)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_diss(6)), &
                     stt=(/nt/), cnt=(/1/))
    !
    do i = 47,49 ! => box-average values of the main COMBINE isotopic variables #1 to #3
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_isot(i-46)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 50,51 ! => box-average values of the main COMBINE isotopic variables #5 to #6
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_isot(i-45)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    ! i=52 corresponds to former unattributed variable #19 (global value)
    !
    i = 53 ! => box-average value of the main COMBINE dissolved variable #8 (SO4^2-)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(vartot_diss(8)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 54
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(ph_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 55
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(temp_tot-273.15), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 56
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(salin_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 57
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(O2_atm_level), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 58
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(CO2_atm_level), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 59
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(O2_atm_conc), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 60
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(CO2_atm_conc), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 61
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(xPOPexport/1.d+12), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 62
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fanthros/1.d+12), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 63
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fco2atm_ocean(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 64
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fco2atm_ocean_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 65
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(finorgC(:)) ,&
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 66
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fdissol_carb(3)*var_part(5,3)*box_vol(3)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 67
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fsilw+fbasw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 68
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fbasw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 69
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fcarbw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 70
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fkerw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 71
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(freef(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 72
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(freef_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 73
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fdep_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 74
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fodc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 75
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fodc_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 76
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fsink_inorg(:)*var_part(5,:)*box_vol(1:nbasin)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 77
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fpw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 78
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fbioC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 79
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(F_seafloor_cdiss(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 80
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(F_seafloor_cdiss_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 81
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(ftrap), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 82
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fCO2_crust(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 83
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSO4_basin(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 84
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSO4_crust(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 85
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(FrivLi), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 86
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(dLiriv), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 87
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(total_cont_POC_export), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 88
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fvol+sum(fmor)+ftrap+fanthros), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 89
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(discharge), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 90
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(tss), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 91
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(epsiC(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 92
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(epsiCont), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 93
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(temp_box(nbasin)), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 94
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(ws(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 95
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(rho_sed*fin_sed(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 96
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(Corg_BE(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 97
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fodp(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 98
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fphos(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 99
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fhydP(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 100
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fcarbsulfw+fsilsulfw+fH2SO4sulfw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 101
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fsilsulfw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 102
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fcarbsulfw), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 103
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fSulfRed(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 104
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fSulfRed_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    i = 105
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fO2_odc(:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    !
    i = 106
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fO2_odc_tot), &
                     stt=(/nt/), cnt=(/1/))
    !
    do i = 107,111 ! => climatic parameter # 1 to 5
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(cpvec(i-106)), &
                     stt=(/nt/), cnt=(/1/))
    end do


  ! Close output file
  ! -----------------

  call close_file(fid)



end subroutine


end module
