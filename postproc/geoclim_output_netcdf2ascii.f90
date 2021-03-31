module local_subroutine
implicit none
contains

  function get_len(string)
    character(len=*), intent(in):: string
    integer:: get_len
    get_len = 0
    do while (string(get_len+1:get_len+1)/=' ')
      get_len = get_len+1
    end do
  end function

  subroutine nf90check(ierr,message)
    use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in), optional:: message
    if (ierr/=NF90_NOERR) then
      if (present(message)) print *, message
      print *, nf90_strerror(ierr)
    end if
  end subroutine

end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program geoclim_output_netcdf2ascii
use netcdf
use local_subroutine, only: get_len, nf90check
implicit none

! ********************************************************************************
  include '/Users/yves/fortran/GEOCLIM4_ber/source/combine_foam.inc'
  character(len=*), parameter:: iofile=geoclim_path//'initfiles/IO_conditions.txt'
  character(len=*), parameter:: outpath='/Users/yves/fortran/GEOCLIM4_ber/OUTPUT/ascii_postprod/' !geoclim_path//'OUTPUT/'
  integer, parameter :: nskip = 99
! ********************************************************************************
  character(len=200) :: dummystr
  character(len=200), dimension(nGEOoutvar):: fname, varname
  character(len=200), dimension(nBIOoutvar):: Bfname, Bvarname
  integer, dimension(nGEOoutvar) :: fnum, fid, varid
  integer, dimension(nBIOoutvar) :: Bfnum, Bfid, Bvarid
  integer :: nl, rnl, opl, ierr, ntime, tfid, tvarid
  double precision:: vector(1)


! ###################
  coupling_biogeo = 0
! ###################


  ! open IO file:
  !--------------

  open(666,status='old',action='read',file=iofile)

  do k = 1,4
    read(666,*)
  end do

  run_name = ' '
  output_path = ' '
  read(666,*) dummystr, run_name
  read(666,*) dummystr, output_path
  rnl = get_len(run_name)
  opl = get_len(output_path)

  do k = 1,nskip-6
    read(666,*)
  end do

    


  ! create ascii output:
  !----------------------

  open (7,file=outpath//'var'//run_name,status='REPLACE')
  open (11,file=outpath//'box1'//run_name,status='REPLACE')
  open (12,file=outpath//'box2'//run_name,status='REPLACE')
  open (13,file=outpath//'box3'//run_name,status='REPLACE')
  open (14,file=outpath//'box4'//run_name,status='REPLACE')
  open (15,file=outpath//'box5'//run_name,status='REPLACE')
  open (16,file=outpath//'box6'//run_name,status='REPLACE')
  open (17,file=outpath//'box7'//run_name,status='REPLACE')
  open (18,file=outpath//'box8'//run_name,status='REPLACE')
  open (19,file=outpath//'box9'//run_name,status='REPLACE')
  open (20,file=outpath//'box10'//run_name,status='REPLACE')
  open (21,file=outpath//'cflux'//run_name,status='REPLACE')
  open (22,file=outpath//'chimie'//run_name,status='REPLACE')
  open (23,file=outpath//'lysocline'//run_name,status='REPLACE')
  open (24,file=outpath//'trap1'//run_name,status='REPLACE')
  open (25,file=outpath//'trap2'//run_name,status='REPLACE')
  open (26,file=outpath//'seacarb_diss'//run_name,status='REPLACE')
  open (27,file=outpath//'speciation'//run_name,status='REPLACE')
  open (28,file=outpath//'speciation_c13'//run_name,status='REPLACE')
  open (29,file=outpath//'forcing'//run_name,status='REPLACE')
  open (205,file=outpath//'lithium'//run_name,status='REPLACE')
  if (coupling_biogeo==1) then
    open (100,file=outpath//'biodiv1'//run_name,status='REPLACE')
    open (101,file=outpath//'biodiv3'//run_name,status='REPLACE')
    open (102,file=outpath//'biodiv6'//run_name,status='REPLACE')
    open (103,file=outpath//'biodiv8'//run_name,status='REPLACE')
    open (204,file=outpath//'species_count'//run_name,status='REPLACE')
  end if



  ! open geoclim netcdf output:
  !----------------------------

  ! get file and variables names:
  do k = 1,nGEOoutvar
    read(666,*) dummystr, fname(k), varname(k)
    print*,k,varname(k)
  end do

  ! loop on all files:

  do i=2,nGEOoutvar
  ! replace "-" by previous line's file name
    dummystr = ' '
    dummystr = fname(i)
    if (dummystr(1:1) == '-') then
      dummystr='#'
      j = i
      do while (dummystr(1:1)=='#')
        j = j-1
        dummystr = fname(j)
      end do
      fname(i) = dummystr
    else
      if (dummystr(1:1) /= '#') then
        ! add output directory path and extension:
        nl = get_len(dummystr)
        fname(i) = geoclim_path//output_path(1:opl)//dummystr(1:nl)//run_name(1:rnl)//'.nc'
      end if
    end if
  end do


  ! file openning:

  fnum(5) = 5 ! the first 4 variables are dimensions variables + volume + surface, and don't need a file
  do k = 5,nGEOoutvar  ! loop on all file names
    i = 1
    do while ( fname(k) /= fname(i) ) ! loop on all previous file names
      i = i + 1
    end do
    fnum(k) = i
  end do
  do k = 1,nGEOoutvar
    dummystr = fname(k)
    if (dummystr(1:1)=='#') then
      ! "#" file name => don't write this variable
      ! => put 0, so that it will not be read
      fnum(k) = 0
    end if
  end do

  k = 0
  do i = 5,nGEOoutvar
    if ( fnum(i) > k ) then ! if "new" file
      ierr = nf90_open( fname(i), NF90_NOWRITE , fid(i) )
      call nf90check(ierr,'openning file '//fname(i))
      ! get time length:
      ierr = nf90_inq_dimid( fid(i) , varname(2), j )
      call nf90check(ierr,'get ID of variable'//varname(2))
      ierr = nf90_inquire_dimension( fid(i) , j , len=ntime )
      call nf90check(ierr)
      tfid = fid(i)
      ierr = nf90_inq_varid( tfid , varname(2) , tvarid )
      call nf90check(ierr)
      k = fnum(i)
    else
      fid(i) = fid(fnum(i))
    end if
  end do

  ! get variables identifier:
  do k = 5,nGEOoutvar
    if (fnum(k) > 0) then
      ierr = nf90_inq_varid( fid(k) , varname(k), varid(k) )
      call nf90check(ierr,'get ID of variable '//varname(k))
    end if
  end do




  ! open biodiv netcdf output:
  !---------------------------

  if (coupling_biogeo==1) then

    do k = 1,3+3+7
      read(666,*)
    end do

    do k = 1,nBIOoutvar
      read(666,*) dummystr, Bfname(k), Bvarname(k)
    end do

    ! loop on all files:

    do i=2,nBIOoutvar
    ! replace "-" by previous line's file name
      dummystr = ' '
      dummystr = Bfname(i)
      if (dummystr(1:1) == '-') then
        dummystr='#'
        j = i
        do while (dummystr(1:1)=='#')
          j = j-1
          dummystr = Bfname(j)
        end do
        Bfname(i) = dummystr
      else
        if (dummystr(1:1) /= '#') then
          ! add output directory path and extension:
          nl = get_len(dummystr)
          Bfname(i) = geoclim_path//output_path(1:opl)//dummystr(1:nl)//run_name(1:rnl)//'.nc'
        end if
      end if
    end do


    ! file openning:

    Bfnum(6) = 6 ! the first 5 variables are dimensions variables + volume, and don't need a file
    do k = 6,nBIOoutvar  ! loop on all file names
      i = 1
      do while ( Bfname(k) /= Bfname(i) ) ! loop on all previous file names
        i = i + 1
      end do
      Bfnum(k) = i
    end do
    do k = 1,nBIOoutvar
      dummystr = Bfname(k)
      if (dummystr(1:1)=='#') then
        ! "#" file name => don't write this variable
        ! => put 0, so that it will not be read
        Bfnum(k) = 0
      end if
    end do

    k = 0
    do i = 6,nBIOoutvar
      if ( Bfnum(i) > k ) then ! if "new" file
        ierr = nf90_open( Bfname(i), NF90_NOWRITE , Bfid(i) )
        call nf90check(ierr,'open file '//Bfname(i))
        k = Bfnum(i)
      else
        Bfid(i) = Bfid(Bfnum(i))
      end if
    end do

    ! get variables identifier:
    do k = 6,nBIOoutvar
      if (Bfnum(k) > 0) then
        ierr = nf90_inq_varid( Bfid(k) , Bvarname(k), Bvarid(k) )
        call nf90check(ierr,'get ID of variable '//Bvarname(k))
      end if
    end do

  end if








  do k = 1,ntime


    ! read geoclim netcdf output:
    !----------------------------

    ierr = nf90_get_var( tfid, tvarid , vector , start=(/k/), count=(/1/) )
    t = vector(1)
    !
    do i = 5,24
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  var(i-4,:)              , start=(/1,k/), count=(/nbasin,1/)  )
    end do
    !
    i = 25
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  h2co3(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 26
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  hco3(:)                 , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 27
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  co3(:)                  , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 28
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  dh2co3(:)               , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 29
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  dhco3(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 30
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  dco3(:)                 , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 31
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  ph(:)                   , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 32
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  omega(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 33
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  temp_box(:)             , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 34
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  salin(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 35
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  dplysc(:)               , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 36
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  dplysa(:)               , start=(/1,k/), count=(/nbasin,1/)  )
    !
    do i = 37,47
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  var_tot(i-36:i-36)       , start=(/k/), count=(/1/)           )
    end do
    !
    do i = 48,50
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  var_tot(i-35:i-35)       , start=(/k/), count=(/1/)           )
    end do
    !
    do i = 51,54
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  var_tot(i-34:i-34)       , start=(/k/), count=(/1/)           )
    end do
    !
    i = 55
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               ph_tot = vector(1)
    !
    i = 56
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               omega_tot = vector(1)
    !
    i = 57
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               temp_tot = vector(1)
    !
    i = 58
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               salin_tot = vector(1)
    !
    i = 59
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               O2_atm_level = vector(1)
    !
    i = 60
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               CO2_atm_level = vector(1)
    !
    i = 63
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               xPOPexport = vector(1)
    !
    i = 64
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fanthros = vector(1)
    !
    i = 65
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fco2atm_ocean(:)        , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 66
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fco2atm_ocean_tot = vector(1)
    !
    i = 67
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  finorgC(:)              , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 68
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fdissol_carb(3:3)       , start=(/k/), count=(/1/)           )
    !
    i = 69
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fsilw = vector(1)
    !
    i = 70
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fbasw = vector(1)
    !
    i = 71
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fcarbw = vector(1)
    !
    i = 72
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fkerw = vector(1)
    !
    i = 73
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  freef(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 74
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               freef_tot = vector(1)
    !
    i = 75
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fodc_tot = vector(1)
    !
    i = 76
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fodc(:)                 , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 77
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fodc_tot = vector(1)
    !
    i = 78
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fsink_inorg(:)          , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 79
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               fpw = vector(1)
    !
    i = 80
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fbioC(:)                , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 81
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  F_seafloor_cdiss(:)     , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 82
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               F_seafloor_cdiss_tot = vector(1)
    !
    i = 83
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               ftrap = vector(1)
    !
    i = 84
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fCO2_crust(:)           , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 85
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fSO4_basin(:)           , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 86
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  fSO4_crust(:)           , start=(/1,k/), count=(/nbasin,1/)  )
    !
    i = 87
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               FrivLi = vector(1)
    !
    i = 88
      vector = 0
      if (fnum(i)>0)  ierr = nf90_get_var(  fid(i), varid(i),  vector                  , start=(/k/), count=(/1/)           )
                                                               dLiriv = vector(1)


    ! read biodiv netcdf output:
    !----------------------------

    if (coupling_biogeo==1) then

      i = 6
      if (Bfnum(i)>0)  ierr = nf90_get_var(  Bfid(i), Bvarid(i),  var_bio     , start=(/1,1,k/), count=(/nequat,nbasin,1/)  )
      !
      i = 7
      if (Bfnum(i)>0)  ierr = nf90_get_var(  Bfid(i), Bvarid(i),  specount             , start=(/1,k/), count=(/nbasin,1/)  )

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
             (fbioC(j),j=1,9)

    write(205,12)t,FrivLi,dLiriv

    write(23,10)t,(dplysc(i),i=1,nbasin-1)

    write(22,15)t,var_tot(1),var_tot(2),ph_tot,omega_tot, &
                    (pH(j),j=1,9),(omega(j),j=1,9)

    write(24,12)t,ftrap,(fCO2_crust(j),j=1,9),(fSO4_basin(j),j=1,9),(fSO4_crust(j),j=1,9)
    write(26,15)t,(F_seafloor_cdiss(j),j=1,nbasin-1)

    write(27,12)t,(h2co3(j),j=1,nbasin-1),(hco3(j),j=1,nbasin-1),(co3(j),j=1,nbasin-1)
    write(28,12)t,(dh2co3(j),j=1,nbasin-1),(dhco3(j),j=1,nbasin-1),(dco3(j),j=1,nbasin-1)

    write(29,17)t/1e6,var(12,nbasin)/0.052d+18,temp_box(1),temp_box(3), &
                  temp_box(6),temp_box(8),salin(3),salin(4),salin(6),dco3(3)*1d3,dco3(6)*1d3

    if (coupling_biogeo==1) then
      write(100,12)(var_bio(i,1),i=1,nequat)
      write(101,12)(var_bio(i,3),i=1,nequat)
      write(102,12)(var_bio(i,6),i=1,nequat)
      write(103,12)(var_bio(i,8),i=1,nequat)
      write(204,19)t/1e6,specount(1),specount(3),specount(6),specount(8)
    end if


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
      ierr = nf90_close( fid(k) )
      k = fnum(i)
    end if
  end do

  if (coupling_biogeo==1) then
    k = 0
    do i = 6,nBIOoutvar
      if ( Bfnum(i) > k ) then ! if "new" file
        ierr = nf90_close( Bfid(k) )
        k = Bfnum(i)
      end if
    end do
  end if



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
  if (coupling_biogeo==1) then
    close(100)
    close(101)
    close(102)
    close(103)
    close(204)
  end if



end program
