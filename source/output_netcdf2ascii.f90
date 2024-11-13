module output_netcdf2ascii_mod
implicit none
contains


subroutine output_netcdf2ascii(COMB_outvar_info)
use netcdf
use netcdf_io_module, only: nf90_check
use io_module, only: netcdf_output_var, DEFAULT_FILLVAL
use constante, only: PI_n_CO2_atm

! ********************************************************************************
  include 'combine_foam.inc'
! ********************************************************************************
  character(len=200) :: dummystr
  integer, dimension(nCOMBoutvar):: varid
  integer:: fid, Bfid
  integer:: ierr, ntime, Bntime, tfid, tvarid, Btvarid
  double precision:: vector(1)
  logical:: open_COMB_file, got_var(nCOMBoutvar)
  ! Declare structures for netCDF output info (cannot be declared in combine_foam.inc)
  type(netcdf_output_var), dimension(nCOMBoutvar), intent(in) :: COMB_outvar_info


  print *
  print *
  print *, '***********************************************'
  print *, '* conversion of netCDF output in ASCII format *'
  print *, '***********************************************'
  print *


  ! create ascii output:
  !----------------------

  open (7,file=trim(output_directory)//'var'//run_name,status='REPLACE')
  open (11,file=trim(output_directory)//'box1'//run_name,status='REPLACE')
  open (12,file=trim(output_directory)//'box2'//run_name,status='REPLACE')
  open (13,file=trim(output_directory)//'box3'//run_name,status='REPLACE')
  open (14,file=trim(output_directory)//'box4'//run_name,status='REPLACE')
  open (15,file=trim(output_directory)//'box5'//run_name,status='REPLACE')
  open (16,file=trim(output_directory)//'box6'//run_name,status='REPLACE')
  open (17,file=trim(output_directory)//'box7'//run_name,status='REPLACE')
  open (18,file=trim(output_directory)//'box8'//run_name,status='REPLACE')
  open (19,file=trim(output_directory)//'box9'//run_name,status='REPLACE')
  open (20,file=trim(output_directory)//'box10'//run_name,status='REPLACE')
  open (21,file=trim(output_directory)//'cflux'//run_name,status='REPLACE')
  open (22,file=trim(output_directory)//'chimie'//run_name,status='REPLACE')
  open (23,file=trim(output_directory)//'lysocline'//run_name,status='REPLACE')
  open (24,file=trim(output_directory)//'trap1'//run_name,status='REPLACE')
  open (25,file=trim(output_directory)//'trap2'//run_name,status='REPLACE')
  open (26,file=trim(output_directory)//'seacarb_diss'//run_name,status='REPLACE')
  open (27,file=trim(output_directory)//'speciation'//run_name,status='REPLACE')
  open (28,file=trim(output_directory)//'speciation_c13'//run_name,status='REPLACE')
  open (29,file=trim(output_directory)//'forcing'//run_name,status='REPLACE')
  open (205,file=trim(output_directory)//'lithium'//run_name,status='REPLACE')



  ! open geoclim netcdf output:
  !----------------------------

  ierr = nf90_open(COMB_ofile_name, NF90_NOWRITE, fid)
  call nf90_check(ierr, 'Error while openning file "'//trim(COMB_ofile_name)// &
                  '". Cannot convert COMBINE netCDF outputs in ASCII format', kill=.false.)

  if (ierr==NF90_NOERR) then
    open_COMB_file = .true.
    ! get time length:
    ierr = nf90_inq_dimid(fid, COMB_time_dimname, tvarid)
    call nf90_check(ierr, 'Error while getting ID of dimension "'//trim(COMB_time_dimname)//'"')
    ierr = nf90_inquire_dimension(fid, tvarid, len=ntime)
    call nf90_check(ierr, 'Error while inquiring length of dimension "'//trim(COMB_time_dimname)//'"')
    ierr = nf90_inq_varid(fid, COMB_time_dimname, tvarid)
    call nf90_check(ierr, 'Error while getting ID of variable "'//trim(COMB_time_dimname)//'"')
    ! get variables identifier:
    do k = 4,nCOMBoutvar ! Skip first 3 variables (box volume, surface and sedimentary surface)
      if (COMB_outvar_info(k)%writevar) then
        ierr = nf90_inq_varid(fid, COMB_outvar_info(k)%vname, varid(k))
        call nf90_check(ierr, 'Error while getting ID of variable "'//trim(COMB_outvar_info(k)%vname)//'"', kill=.false.)
        got_var(k) = .true.
      else
        got_var(k) = .false. ! indicate that variable was not saved
      end if
    end do
  else
    open_COMB_file = .false. ! indicate that netCDF file couldn't be open
  end if



  ! ========================================================================== !

  ! Fill-value on all variables (in case they were not saved)
  t                    = DEFAULT_FILLVAL
  var                  = DEFAULT_FILLVAL
  h2co3(:)             = DEFAULT_FILLVAL
  hco3(:)              = DEFAULT_FILLVAL
  co3(:)               = DEFAULT_FILLVAL
  dh2co3(:)            = DEFAULT_FILLVAL
  dhco3(:)             = DEFAULT_FILLVAL
  dco3(:)              = DEFAULT_FILLVAL
  ph(:)                = DEFAULT_FILLVAL
  omega(:)             = DEFAULT_FILLVAL
  temp_box(:)          = DEFAULT_FILLVAL
  salin(:)             = DEFAULT_FILLVAL
  dplysc(:)            = DEFAULT_FILLVAL
  dplysa(:)            = DEFAULT_FILLVAL
  var_tot              = DEFAULT_FILLVAL
  ph_tot               = DEFAULT_FILLVAL
  temp_tot             = DEFAULT_FILLVAL
  salin_tot            = DEFAULT_FILLVAL
  O2_atm_level         = DEFAULT_FILLVAL
  CO2_atm_level        = DEFAULT_FILLVAL
  xPOPexport           = DEFAULT_FILLVAL
  fanthros             = DEFAULT_FILLVAL
  fco2atm_ocean(:)     = DEFAULT_FILLVAL
  fco2atm_ocean_tot    = DEFAULT_FILLVAL
  finorgC(:)           = DEFAULT_FILLVAL
  fdissol_carb(:)      = DEFAULT_FILLVAL
  fsilw                = DEFAULT_FILLVAL
  fbasw                = DEFAULT_FILLVAL
  fkerw                = DEFAULT_FILLVAL
  freef(:)             = DEFAULT_FILLVAL
  freef_tot            = DEFAULT_FILLVAL
  fdep_tot             = DEFAULT_FILLVAL
  fodc(:)              = DEFAULT_FILLVAL
  fodc_tot             = DEFAULT_FILLVAL
  fsink_inorg(:)       = DEFAULT_FILLVAL
  fpw                  = DEFAULT_FILLVAL
  fbioC(:)             = DEFAULT_FILLVAL
  F_seafloor_cdiss(:)  = DEFAULT_FILLVAL
  F_seafloor_cdiss_tot = DEFAULT_FILLVAL
  ftrap                = DEFAULT_FILLVAL
  fCO2_crust(:)        = DEFAULT_FILLVAL
  fSO4_basin(:)        = DEFAULT_FILLVAL
  fSO4_crust(:)        = DEFAULT_FILLVAL
  FrivLi               = DEFAULT_FILLVAL
  dLiriv               = DEFAULT_FILLVAL


  do k = 1,ntime

    ! read geoclim netcdf output:
    !----------------------------

    if (open_COMB_file) then
      ierr = nf90_get_var( fid, tvarid , vector , start=(/k/), count=(/1/) )
      t = vector(1)
      !
      do i = 4,23
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var(i-3,:)           , start=(/1,k/), count=(/nbasin,1/)  )
      end do
      !
      i = 24
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  h2co3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 25
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  hco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 26
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  co3(:)               , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 27
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dh2co3(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 28
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dhco3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 29
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 30
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  ph(:)                , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 31
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  omega(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 32
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  temp_box(:)          , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 33
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  salin(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 34
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dplysc(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 35
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  dplysa(:)            , start=(/1,k/), count=(/nbasin,1/)  )
      !
      do i = 36,46
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_tot(i-35:i-35)    , start=(/k/), count=(/1/)          )
      end do
      !
      do i = 47,49
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_tot(i-34:i-34)    , start=(/k/), count=(/1/)          )
      end do
      !
      do i = 50,53
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  var_tot(i-33:i-33)    , start=(/k/), count=(/1/)          )
      end do
      !
      i = 54
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         ph_tot = vector(1)
      !
      i = 55
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         temp_tot = vector(1)
      !
      i = 56
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         salin_tot = vector(1)
      !
      i = 57
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         O2_atm_level = vector(1)
      !
      i = 58
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         CO2_atm_level = vector(1)
      !
      i = 61
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         xPOPexport = vector(1)
      !
      i = 62
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fanthros = vector(1)
      !
      i = 63
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fco2atm_ocean(:)     , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 64
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fco2atm_ocean_tot = vector(1)
      !
      i = 65
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  finorgC(:)           , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 66
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fdissol_carb(3:3)    , start=(/k/), count=(/1/)           )
      !
      i = 67
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fsilw = vector(1)
      !
      i = 68
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fbasw = vector(1)
      !
      i = 69
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fcarbw = vector(1)
      !
      i = 70
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fkerw = vector(1)
      !
      i = 71
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  freef(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 72
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         freef_tot = vector(1)
      !
      i = 73
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fdep_tot = vector(1)
      !
      i = 74
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fodc(:)              , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 75
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fodc_tot = vector(1)
      !
      i = 76
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fsink_inorg(:)       , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 77
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         fpw = vector(1)
      !
      i = 78
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fbioC(:)             , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 79
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  F_seafloor_cdiss(:)  , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 80
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         F_seafloor_cdiss_tot = vector(1)
      !
      i = 81
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         ftrap = vector(1)
      !
      i = 82
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fCO2_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 83
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fSO4_basin(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 84
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  fSO4_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
      !
      i = 85
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         FrivLi = vector(1)
      !
      i = 86
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         dLiriv = vector(1)
      i = 87
        vector = 0
        if (got_var(i))  ierr = nf90_get_var(  fid, varid(i),  vector               , start=(/k/), count=(/1/)           )
                         total_cont_POC_export = vector(1)
    end if



    ! write ascii output:
    !--------------------

    write(7,9)t,O2_atm_level,CO2_atm_level, &
                  fanthros,fco2atm_ocean(1),fco2atm_ocean(3), &
                  fco2atm_ocean(6),fco2atm_ocean(8),xPOPexport, &
                  finorgC(3),fdissol_carb(3),dplysa(3)

    write(11,11)t,(var(j,1),j=1,nvar_real)
    write(12,11)t,(var(j,2),j=1,nvar_real)
    write(13,11)t,(var(j,3),j=1,nvar_real)
    write(14,11)t,(var(j,4),j=1,nvar_real)
    write(15,11)t,(var(j,5),j=1,nvar_real)
    write(16,11)t,(var(j,6),j=1,nvar_real)
    write(17,11)t,(var(j,7),j=1,nvar_real)
    write(18,11)t,(var(j,8),j=1,nvar_real)
    write(19,11)t,(var(j,9),j=1,nvar_real)
    write(20,11)t,(var(j,10),j=1,nvar_real)

    write(21,12)t,fsilw,fbasw,fcarbw,fkerw,freef_tot, &
                    fdep_tot,fodc_tot, &
                    (freef(j),j=1,9), &
                    (fodc(j),j=1,9), &
             (fsink_inorg(j)*var(10,j),j=1,9),fpw, &
             (fbioC(j),j=1,9),total_cont_POC_export

    write(205,12)t,FrivLi,dLiriv

    write(23,10)t,(dplysc(i),i=1,nbasin-1)

    write(22,15)t,var_tot(1),var_tot(2),ph_tot, &
                    (pH(j),j=1,9),(omega(j),j=1,9)

    write(24,12)t,ftrap,(fCO2_crust(j),j=1,9),(fSO4_basin(j),j=1,9),(fSO4_crust(j),j=1,9)
    write(26,15)t,(F_seafloor_cdiss(j),j=1,nbasin-1)

    write(27,12)t,(h2co3(j),j=1,nbasin-1),(hco3(j),j=1,nbasin-1),(co3(j),j=1,nbasin-1)
    write(28,12)t,(dh2co3(j),j=1,nbasin-1),(dhco3(j),j=1,nbasin-1),(dco3(j),j=1,nbasin-1)

    write(29,17)t/1e6,var(12,nbasin)/PI_n_CO2_atm,temp_box(1),temp_box(3), &
                  temp_box(6),temp_box(8),salin(3),salin(4),salin(6),dco3(3)*1d3,dco3(6)*1d3


   9    format(1f15.5,1x,2(f15.6,1x),20(es15.6e3,1x))
   11   format(21(es15.6e3,1x))
   12   format(500(es15.6e3,1x))
   10   format(10(es15.6e3,1x))
   15   format(23(es15.6e3,1x))
   16   format(1921(es15.6e3,1x))
   17   format(30(f10.4,1x))
   18   format(1920(es15.6e2,1x))
   19   format(1f15.5,1x,4i5)



  end do



  ! close netcdf output:
  !---------------------

  if (open_COMB_file) ierr = nf90_close(fid)



  ! close ascii output:
  !--------------------

  close(7)
  close(11)
  close(12)
  close(13)
  close(14)
  close(15)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)
  close(23)
  close(24)
  close(25)
  close(26)
  close(27)
  close(28)
  close(29)
  close(205)



end subroutine


end module
