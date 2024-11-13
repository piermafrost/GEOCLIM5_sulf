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
        xPOPexport=xPOPexport + fsink(j)*var(7,j)*vol(j)
        !totalCO2_exch=totalCO2_exch+fco2atm_ocean(j)
    enddo

    vol_tot=0.
    var_tot(:)=0.
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
        vol_tot              = vol_tot + vol(j)
        var_tot(1:11)        = var_tot(1:11) + var(1:11,j)*vol(j)              !| concentration
        var_tot(20)          = var_tot(20) + var(20,j)*vol(j)                  !| of elements
        var_tot(13)          = var_tot(13) + var(13,j)*var(1,j)*vol(j)         ! DIC d13C
        var_tot(14)          = var_tot(14) + var(14,j)*var(10,j)*vol(j)        ! PIC d13C
        var_tot(15)          = var_tot(15) + var(15,j)*var(9,j)*vol(j)         ! POC d13C
        var_tot(17:18)       = var_tot(17:18) + var(17:18,j)*var(5:6,j)*vol(j) ! diss & PIC Sr iso ration
        var_tot(19)          = var_tot(19) + var(19,j)*var(20,j)*vol(j)        ! d7Li
        ph_tot               = ph_tot + pH(j)*vol(j)
        temp_tot             = temp_tot + temp_box(j)*vol(j)
        salin_tot            = salin_tot + salin(j)*vol(j)
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
        fdep_tot = fdep_tot + fsink_inorg(j)*var(10,j)*(vol(j))
    end do
    do j0=1,nappcont-1!<- skip atmospheric box, counted in "appcont"
        j = jbox_appcont(j0)
        fdep_tot = fdep_tot - (1-gotoshelf)*fsink_inorg(j)*var(10,j)*(vol(j))
    end do

    ! Normalize mean isotopic ratios
    var_tot(13)    = var_tot(13) / var_tot(1)      ! DIC d13C
    !
    if (var_tot(10) == 0d0) then
        var_tot(14) = COMB_outvar_info(14)%fillval
    else
        var_tot(14) = var_tot(14) / var_tot(10)    ! PIC d13C
    end if
    !
    var_tot(15)    = var_tot(15) / var_tot(9)      ! POC d13C
    !
    if (var_tot(5) == 0d0) then
        var_tot(17) = COMB_outvar_info(17)%fillval
    else
        var_tot(17) = var_tot(17) / var_tot(5)     ! diss Sr iso ration
    end if
    !
    if (var_tot(6) == 0d0) then
        var_tot(18) = COMB_outvar_info(18)%fillval
    else
        var_tot(18) = var_tot(18) / var_tot(6)     ! PIC Sr iso ration
    end if
    !
    var_tot(19)    = var_tot(19) / var_tot(20)     ! d7Li

    ! Normalize mean intensive variables (e.g., concentration)
    var_tot(1:11)  = var_tot(1:11) / vol_tot       !| concentration
    var_tot(20)    = var_tot(20) / vol_tot         !| of elements
    ph_tot         = ph_tot / vol_tot
    temp_tot       = temp_tot / vol_tot
    salin_tot      = salin_tot / vol_tot

    ! Present-day normalized mol of O2 and CO2
    O2_atm_level  = var(11,nbasin) / PI_n_O2_atm
    CO2_atm_level = var(12,nbasin) / PI_n_CO2_atm

    ! convert mol to volume fraction (= molar fraction) for 02 and CO2
    O2_atm_conc   =  ( var(11,nbasin) / ( PI_n_rest_of_atm + var(11,nbasin) + var(12,nbasin) ) )*100 ! (%)
    CO2_atm_conc  =  ( var(12,nbasin) / ( PI_n_rest_of_atm + var(11,nbasin) + var(12,nbasin) ) )*1e6 ! (ppm)



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

    do i = 4,23 ! => the main COMBINE variables #1 to #20 (see creades.f)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(var(i-3,:)), &
                     stt=(/1,nt/), cnt=(/nbasin,1/))
    end do
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
    do i = 36,46 ! => box-average values of the main COMBINE variables #1 to #11 (see creades.f)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(var_tot(i-35)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 47,49 ! => box-average values of the main COMBINE variables #13 to #15 (see creades.f)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(var_tot(i-34)), &
                     stt=(/nt/), cnt=(/1/))
    end do
    !
    do i = 50,53 ! => box-average values of the main COMBINE variables #17 to #20 (see creades.f)
      if (COMB_outvar_info(i)%writevar) &
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(var_tot(i-33)), &
                     stt=(/nt/), cnt=(/1/))
    end do
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
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real0D=real(fdissol_carb(3)*var(10,3)*vol(3)), &
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
        call put_var(fid, varname=COMB_outvar_info(i)%vname, var_real1D=real(fsink_inorg(:)*var(10,:)*vol(1:nbasin)), &
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
