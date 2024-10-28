module output_netcdf2ascii_mod
implicit none
contains

  subroutine nf90check(ierr,message)
    use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in), optional:: message
    if (ierr/=NF90_NOERR) then
      if (present(message)) print *, message
      print *, nf90_strerror(ierr)
    end if
  end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine output_netcdf2ascii()
use netcdf
use constante, only: PI_n_CO2_atm

! ********************************************************************************
  include 'combine_foam.inc'
! ********************************************************************************
  character(len=200) :: dummystr
  integer, dimension(nGEOoutvar):: nt, varid
  integer:: fnum(nGEOoutvar)
  integer :: nl, rnl, ierr, ntime, tfid, tvarid
  double precision:: vector(1)


  ! create ascii output:
  !----------------------

  open (7,file=trim(output_path)//'var'//run_name,status='REPLACE')
  open (11,file=trim(output_path)//'box1'//run_name,status='REPLACE')
  open (12,file=trim(output_path)//'box2'//run_name,status='REPLACE')
  open (13,file=trim(output_path)//'box3'//run_name,status='REPLACE')
  open (14,file=trim(output_path)//'box4'//run_name,status='REPLACE')
  open (15,file=trim(output_path)//'box5'//run_name,status='REPLACE')
  open (16,file=trim(output_path)//'box6'//run_name,status='REPLACE')
  open (17,file=trim(output_path)//'box7'//run_name,status='REPLACE')
  open (18,file=trim(output_path)//'box8'//run_name,status='REPLACE')
  open (19,file=trim(output_path)//'box9'//run_name,status='REPLACE')
  open (20,file=trim(output_path)//'box10'//run_name,status='REPLACE')
  open (21,file=trim(output_path)//'cflux'//run_name,status='REPLACE')
  open (22,file=trim(output_path)//'chimie'//run_name,status='REPLACE')
  open (23,file=trim(output_path)//'lysocline'//run_name,status='REPLACE')
  open (24,file=trim(output_path)//'trap1'//run_name,status='REPLACE')
  open (25,file=trim(output_path)//'trap2'//run_name,status='REPLACE')
  open (26,file=trim(output_path)//'seacarb_diss'//run_name,status='REPLACE')
  open (27,file=trim(output_path)//'speciation'//run_name,status='REPLACE')
  open (28,file=trim(output_path)//'speciation_c13'//run_name,status='REPLACE')
  open (29,file=trim(output_path)//'forcing'//run_name,status='REPLACE')
  open (205,file=trim(output_path)//'lithium'//run_name,status='REPLACE')



  ! open geoclim netcdf output:
  !----------------------------

  fnum = GEO_ofile_num

  k = 0
  do i = 5,nGEOoutvar
    if ( fnum(i) > k ) then ! if "new" file
      ierr = nf90_open( GEO_ofile_name(i), NF90_NOWRITE , GEO_ofile_num(i) )
      call nf90check(ierr,'openning file '//GEO_ofile_name(i))
      ! get time length:
      ierr = nf90_inq_dimid( GEO_ofile_num(i) , GEO_varout_name(2), j )
      call nf90check(ierr,'get ID of variable'//GEO_varout_name(2))
      ierr = nf90_inquire_dimension( GEO_ofile_num(i) , j , len=ntime )
      call nf90check(ierr)
      tfid = GEO_ofile_num(i)
      ierr = nf90_inq_varid( tfid , GEO_varout_name(2) , tvarid )
      call nf90check(ierr)
      k = fnum(i)
    else
      GEO_ofile_num(i) = GEO_ofile_num(fnum(i))
    end if
  end do

  ! get variables identifier:
  do k = 5,nGEOoutvar
    if (fnum(k) > 0) then
      ierr = nf90_inq_varid( GEO_ofile_num(k) , GEO_varout_name(k), varid(k) )
      call nf90check(ierr,'get ID of variable '//GEO_varout_name(k))
    end if
  end do



  do k = 1,ntime


    ! read geoclim netcdf output:
    !----------------------------

    ierr = nf90_get_var( tfid, tvarid , vector , start=(/k/), count=(/1/) )
    t = vector(1)
    !
    do i = 6,25
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  var(i-5,:)           , start=(/1,k/), count=(/nbasin,1/)  )
    end do
    !
    i = 26
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  h2co3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 27
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  hco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 28
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  co3(:)               , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 29
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  dh2co3(:)            , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 30
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  dhco3(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 31
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  dco3(:)              , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 32
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  ph(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 33
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  omega(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 34
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  temp_box(:)          , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 35
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  salin(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 36
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  dplysc(:)            , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 37
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  dplysa(:)            , start=(/1,k/), count=(/nbasin,1/)  )
    !
    do i = 38,48
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  var_tot(i-37:i-37)    , start=(/k/), count=(/1/)          )
    end do
    !
    do i = 49,51
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  var_tot(i-36:i-36)    , start=(/k/), count=(/1/)          )
    end do
    !
    do i = 52,55
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  var_tot(i-35:i-35)    , start=(/k/), count=(/1/)          )
    end do
    !
    i = 56
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               ph_tot = vector(1)
    !
    i = 57
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               omega_tot = vector(1)
    !
    i = 58
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               temp_tot = vector(1)
    !
    i = 59
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               salin_tot = vector(1)
    !
    i = 60
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               O2_atm_level = vector(1)
    !
    i = 61
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               CO2_atm_level = vector(1)
    !
    i = 64
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               xPOPexport = vector(1)
    !
    i = 65
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fanthros = vector(1)
    !
    i = 66
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fco2atm_ocean(:)     , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 67
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fco2atm_ocean_tot = vector(1)
    !
    i = 68
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  finorgC(:)           , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 69
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fdissol_carb(3:3)    , start=(/k/), count=(/1/)           )
    !
    i = 70
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fsilw = vector(1)
    !
    i = 71
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fbasw = vector(1)
    !
    i = 72
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fcarbw = vector(1)
    !
    i = 73
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fkerw = vector(1)
    !
    i = 74
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  freef(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 75
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               freef_tot = vector(1)
    !
    i = 76
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fodc_tot = vector(1)
    !
    i = 77
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fodc(:)              , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 78
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fodc_tot = vector(1)
    !
    i = 79
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fsink_inorg(:)       , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 80
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               fpw = vector(1)
    !
    i = 81
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fbioC(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 82
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  F_seafloor_cdiss(:)  , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 83
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               F_seafloor_cdiss_tot = vector(1)
    !
    i = 84
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               ftrap = vector(1)
    !
    i = 85
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fCO2_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 86
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fSO4_basin(:)        , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 87
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  fSO4_crust(:)        , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 88
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               FrivLi = vector(1)
    !
    i = 89
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               dLiriv = vector(1)
    i = 90
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  GEO_ofile_num(i), varid(i),  vector               , start=(/k/), count=(/1/)           )
                                                               total_cont_POC_export = vector(1)



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

    write(22,15)t,var_tot(1),var_tot(2),ph_tot,omega_tot, &
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

  k = 0
  do i = 5,nGEOoutvar
    if ( fnum(i) > k ) then ! if "new" file
      ierr = nf90_close( GEO_ofile_num(k) )
      k = fnum(i)
    end if
  end do



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
