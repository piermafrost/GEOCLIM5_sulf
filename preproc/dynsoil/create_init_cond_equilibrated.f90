module isitiok_mod
implicit none
contains
  subroutine isitok(ierr,message)
    use netcdf
    integer, intent(in):: ierr
    character(len=*), intent(in):: message
    if (ierr/=0) then
      print *, message
      print *, nf90_strerror(ierr)
      stop
    end if
  end subroutine
end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module get_len_mod
implicit none
contains
  function get_len(string)
    character(len=*), intent(in):: string
    integer:: get_len
    integer::k
    get_len = len(string)
    k = 0
    do while (k<get_len)
      if (string(k+1:k+1)/=' ') then
        k = k+1
      else
        get_len = k
      end if
    end do
  end function
end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module find_co2_interval_mod
implicit none
contains
  subroutine find_co2_interval( co2lev, kinf, ksup, xi )
    double precision, intent(in)::  co2lev
    integer, intent(out)::          kinf, ksup
    double precision, intent(out):: xi
    include 'combine_foam.inc'
    kinf = 1
    ksup = nclimber
    ! find by dichotomy:
    do while ( ksup-kinf > 1 )
      k = kinf + (ksup-kinf)/2
      if (co2lev >= co2climber(k)) then
        kinf = k
      else
        ksup = k
      end if
    end do
    ! interpolation coefficient:
    xi = ( co2lev - co2climber(kinf) )  /  ( co2climber(ksup) - co2climber(kinf) )
  end subroutine
end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module dynsoil_interface
implicit none
contains

  function DS_sigma(j)
    integer, intent(in):: j
    double precision:: DS_sigma
    include 'dynsoil_physical_parameters.inc'
    DS_sigma = sigma(j)
  end function

  function DS_scaling_factor()
    double precision:: DS_scaling_factor
    include 'dynsoil_physical_parameters.inc'
    DS_scaling_factor = scaling_factor
  end function

  function DS_epsl()
    double precision:: DS_epsl
    include 'dynsoil_physical_parameters.inc'
    DS_epsl = epsl
  end function

end module


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program create_initial_condition_equilibrated

use netcdf
use isitiok_mod, only: isitok
use get_len_mod, only: get_len
use find_co2_interval_mod, only: find_co2_interval
use dynsoil_interface, only: DS_sigma, DS_scaling_factor, DS_epsl
use dynsoil_empirical_laws, only: erosion, reg_prod_opt, dissolution_constant, eq_reg_thick!, soil_prod_function
implicit none


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  character(len=*), parameter:: ofname='dynsoil_restart_equilib.nc'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  character(len=200):: areafile, lithfile, tempfile, rnfffile, slopfile
  double precision:: CO2_level

  real, parameter:: fillvalue = 9.96921e+36
  include 'combine_foam.inc'
  character(len=200):: area_fname
  double precision:: xi, loc_minslope, Kmain, loc_Esum, loc_Wsum, loc_vol_Wsum
  integer:: fileid, varid, ierr, Odimids(4), Ovarids(9)
  logical:: got_grid


!=========================================================================!
! Argument scanning
!=========================================================================!

  print *
  print *
  print *, 'All input find are assumed to be in ascii format unless their extension is .nc'
  print *, 'Exception: the lithology file must be at netCDF format.'
  print *, 'Temperature and runoff are linearly interpolated to get their values ar the'
  print *, 'given CO2 level (whose units should match the one of climate input files)'
  print *, 'Enter "0" if the climate fields are at a unique CO2 level'
  print *
  print *


  n = command_argument_count()

  write(unit=*, fmt='(A24)', advance='no') 'area input file:        '
  if (n>=1) then
    call get_command_argument(1, areafile)
    print *, areafile(1:get_len(areafile))
  else
    read(unit=*,fmt='(A)') areafile
  end if 

  write(unit=*, fmt='(A24)', advance='no') 'lithology input file:   '
  if (n>=2) then
    call get_command_argument(2, lithfile)
    print *, lithfile(1:get_len(lithfile))
  else
    read(unit=*,fmt='(A)') lithfile
  end if 

  write(unit=*, fmt='(A24)', advance='no') 'temperature input file: '
  if (n>=3) then
    call get_command_argument(3, tempfile)
    print *, tempfile(1:get_len(tempfile))
  else
    read(unit=*,fmt='(A)') tempfile
  end if 

  write(unit=*, fmt='(A24)', advance='no') 'runoff input file:      '
  if (n>=4) then
    call get_command_argument(4, rnfffile)
    print *, rnfffile(1:get_len(rnfffile))
  else
    read(unit=*,fmt='(A)') rnfffile
  end if 

  write(unit=*, fmt='(A24)', advance='no') 'slope input file:       '
  if (n>=5) then
    call get_command_argument(5, slopfile)
    print *, slopfile(1:get_len(slopfile))
  else
    read(unit=*,fmt='(A)') slopfile
  end if 

  write(unit=*, fmt='(A11)', advance='no') 'CO2 level: '
  if (n>=6) then
    call get_command_argument(6, dummychar)
    read(dummychar,*) CO2_level
    print *, CO2_level
  else
    read *, CO2_level
  end if 

  print *
  print *



!=========================================================================!
! read input
!=========================================================================!

got_grid = .false.

! area:
  n = get_len(areafile)

  if (areafile(n-2:n)=='.nc') then ! netCDF format


  ! File openning:
    ierr = nf90_open(areafile, NF90_NOWRITE, fileid)
    call isitok(ierr,'Error while openning input area file '//areafile(1:n))

  ! Loading variable
    ierr = nf90_inq_varid( fileid, 'area', varid )
    call isitok(ierr,'Error while getting "area" variable identifiant in file'//areafile(1:))
    do i=1,nlat
      ierr = nf90_get_var( fileid, varid, areaclimber( 1+nlon*(i-1) : nlon*i ), (/1,i/), (/nlon,1/) )
      call isitok(ierr,'Error while getting variable "area" in file '//areafile(1:n))
    end do

  ! Loading lon and lat variable
    ierr = nf90_inq_varid( fileid, 'lon', varid )
    call isitok(ierr,'Error while getting "lon" variable identifiant in file'//areafile(1:n))
    ierr = nf90_get_var( fileid, varid, ref_x_axis )
    call isitok(ierr,'Error while getting variable "lon" in file'//areafile(1:n))
    ierr = nf90_inq_varid( fileid, 'lat', varid )
    call isitok(ierr,'Error while getting "lat" variable identifiant in file'//areafile(1:n))
    ierr = nf90_get_var( fileid, varid, ref_y_axis )
    call isitok(ierr,'Error while getting variable "lat" in file'//areafile(1:n))
    got_grid = .true.

  ! File closing:
    ierr = nf90_close(fileid)


  else ! "Yves" ascii format assumed

    open( unit=1, file=areafile(1:n), status='old', action='read' )

    do i = 1,npixel
      read(unit=1,fmt=*) areaclimber(i)
    end do
    areaclimber = 1.d12*areaclimber ! in "Yves" ascii format, area is in 1e6km2 instead of m2

    close(unit=1)

  end if

!-------------------------------------------------------------------------!

! lithological fraction:

n = get_len(lithfile)

! File openning:
  ierr = nf90_open(lithfile(1:n), NF90_NOWRITE, fileid)
  call isitok(ierr,'Error while openning input lithology file '//lithfile(1:n))

! Loading variable
  ierr = nf90_inq_varid( fileid, 'litho_frac', varid )
  call isitok(ierr,'Error while getting "litho_frac" variable identifiant in file '//lithfile(1:n))
  do i=1,nlat
    do j=1,nlitho
      ierr = nf90_get_var( fileid, varid, litho_frac( j , 1+nlon*(i-1) : nlon*i ), (/1,i,j/), (/nlon,1,1/) )
      call isitok(ierr,'Error while getting variable "litho_frac" in file '//lithfile(1:n))
    end do
  end do

! File closing:
  ierr = nf90_close(fileid)

!-------------------------------------------------------------------------!

! Climate fields:
!----------------

! temperature: (degrees celsius)
  n = get_len(tempfile)

  if (tempfile(n-2:n)=='.nc') then ! netCDF format


! File openning:
  ierr = nf90_open(tempfile, NF90_NOWRITE, fileid)
  call isitok(ierr,'Error while openning input temp file '//tempfile(1:n))


    if (CO2_level>0) then


    ! Loading variable
      ierr = nf90_inq_varid( fileid, 'CO2', varid )
      call isitok(ierr,'Error while getting "CO2" variable identifiant in file '//tempfile(1:n))
      ierr = nf90_get_var( fileid, varid, co2climber )
      call isitok(ierr,'Error while getting variable "CO2" in file'//tempfile(1:n))
      ierr = nf90_inq_varid( fileid, 'temperature', varid )
      call isitok(ierr,'Error while getting "temperature" variable identifiant in file'//tempfile(1:n))
      do i=1,nlat
        ierr = nf90_get_var( fileid, varid, Tairclimber( 1+nlon*(i-1) : nlon*i , : ), (/1,i,1/), (/nlon,1,nclimber/) )
        call isitok(ierr,'Error while getting variable "temperature" in file'//tempfile(1:n))
      end do


    else


    ! Loading variable
      ierr = nf90_inq_varid( fileid, 'temperature', varid )
      call isitok(ierr,'Error while getting "temperature" variable identifiant in file '//tempfile(1:n))
      do i=1,nlat
        ierr = nf90_get_var( fileid, varid, Tclim( 1+nlon*(i-1) : nlon*i ), (/1,i/), (/nlon,1/) )
        call isitok(ierr,'Error while getting variable "temperature" in file '//tempfile(1:n))
      end do


    end if


    if (.not. got_grid) then
    ! Loading lon and lat variable
      ierr = nf90_inq_varid( fileid, 'lon', varid )
      call isitok(ierr,'Error while getting "lon" variable identifiant in file'//tempfile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_x_axis )
      call isitok(ierr,'Error while getting variable "lon" in file'//tempfile(1:n))
      ierr = nf90_inq_varid( fileid, 'lat', varid )
      call isitok(ierr,'Error while getting "lat" variable identifiant in file'//tempfile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_y_axis )
      call isitok(ierr,'Error while getting variable "lat" in file'//tempfile(1:n))
      got_grid = .true.
    end if

  ! File closing:
    ierr = nf90_close(fileid)



  else ! "Yves" ascii format assumed

    open( unit=1, file=tempfile(1:n), status='old', action='read' )

    if (CO2_level>0) then

      do j=1,nclimber
        read(1,*)co2climber(nclimber+1-j), (Tairclimber(i,nclimber+1-j),i=1,npixel)
      end do

    else

      do i=1,npixel
        read(unit=1,fmt=*) Tclim(i)
      end do

    end if

    close(unit=1)

  end if

!-------------------------------------------------------------------------!

! runoff: (cm/y)
  n = get_len(rnfffile)

  if (rnfffile(n-2:n)=='.nc') then ! netCDF format


  ! File openning:
    ierr = nf90_open(rnfffile, NF90_NOWRITE, fileid)
    call isitok(ierr,'Error while openning input runoff file '//rnfffile(1:n))

    if (CO2_level>0) then

    ! Loading variable
      ierr = nf90_inq_varid( fileid, 'runoff', varid )
      call isitok(ierr,'Error while getting "runoff" variable identifiant in file '//rnfffile(1:n))
      do i=1,nlat
        ierr = nf90_get_var( fileid, varid, Runclimber( 1+nlon*(i-1) : nlon*i , : ), (/1,i,1/), (/nlon,1,nclimber/) )
        call isitok(ierr,'Error while getting variable "runoff" in file'//rnfffile(1:n))
      end do
      Runclimber = 100*Runclimber ! netCDF file should have runoff in m/y


    else

    ! Loading variable
      ierr = nf90_inq_varid( fileid, 'runoff', varid )
      call isitok(ierr,'Error while getting "runoff" variable identifiant in file '//rnfffile(1:n))
      do i=1,nlat
        ierr = nf90_get_var( fileid, varid, runclim( 1+nlon*(i-1) : nlon*i ), (/1,i/), (/nlon,1/) )
        call isitok(ierr,'Error while getting variable "runoff" in file '//rnfffile(1:n))
      end do
      runclim = 100*Runclim ! netCDF file should have runoff in m/y

    end if

    if (.not. got_grid) then
    ! Loading lon and lat variable
      ierr = nf90_inq_varid( fileid, 'lon', varid )
      call isitok(ierr,'Error while getting "lon" variable identifiant in file'//rnfffile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_x_axis )
      call isitok(ierr,'Error while getting variable "lon" in file'//rnfffile(1:n))
      ierr = nf90_inq_varid( fileid, 'lat', varid )
      call isitok(ierr,'Error while getting "lat" variable identifiant in file'//rnfffile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_y_axis )
      call isitok(ierr,'Error while getting variable "lat" in file'//rnfffile(1:n))
      got_grid = .true.
    end if

  ! File closing:
    ierr = nf90_close(fileid)



  else ! "Yves" ascii format assumed

    open( unit=1, file=rnfffile(1:n), status='old', action='read' )

    if (CO2_level>0) then

      do j=1,nclimber
        read(1,*)dummy,(Runclimber(i,nclimber+1-j),i=1,npixel)
      end do

    else

      do i=1,npixel
        read(unit=1,fmt=*) runclim(i)
      end do

    end if

    close(unit=1)

  end if

!-------------------------------------------------------------------------!

! slope:
  n = get_len(slopfile)

  if (slopfile(n-2:n)=='.nc') then ! netCDF format


  ! File openning:
    ierr = nf90_open(slopfile, NF90_NOWRITE, fileid)
    call isitok(ierr,'Error while openning input slope file ')

  ! Loading variable
    ierr = nf90_inq_varid( fileid, 'slope', varid )
    call isitok(ierr,'Error while getting "slope" variable identifiant in file ')
    do i=1,nlat
      ierr = nf90_get_var( fileid, varid, slope( 1+nlon*(i-1) : nlon*i ), (/1,i/), (/nlon,1/) )
      call isitok(ierr,'Error while getting variable "slope" in file ')
    end do

    if (.not. got_grid) then
    ! Loading lon and lat variable
      ierr = nf90_inq_varid( fileid, 'lon', varid )
      call isitok(ierr,'Error while getting "lon" variable identifiant in file'//slopfile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_x_axis )
      call isitok(ierr,'Error while getting variable "lon" in file'//slopfile(1:n))
      ierr = nf90_inq_varid( fileid, 'lat', varid )
      call isitok(ierr,'Error while getting "lat" variable identifiant in file'//slopfile(1:n))
      ierr = nf90_get_var( fileid, varid, ref_y_axis )
      call isitok(ierr,'Error while getting variable "lat" in file'//slopfile(1:n))
      got_grid = .true.
    end if

  ! File closing:
    ierr = nf90_close(fileid)


  else ! "Yves" ascii format assumed

    open( unit=1, file=slopfile(1:n), status='old', action='read' )

    do i = 1,npixel
      read(unit=1,fmt=*) slope(i)
    end do

    close(unit=1)

  end if

!-------------------------------------------------------------------------!

! longitude and latitude

  if (.not. got_grid) then

    n = get_len(lithfile)

  ! File openning:
    ierr = nf90_open(lithfile(1:n), NF90_NOWRITE, fileid)
    call isitok(ierr,'Error while openning input file '//lithfile(1:n))

  ! Loading variable
    ierr = nf90_inq_varid( fileid, 'lon', varid )
    call isitok(ierr,'Error while getting "lon" variable identifiant in file '//lithfile(1:n))
    ierr = nf90_get_var( fileid, varid, ref_x_axis )
    call isitok(ierr,'Error while getting variable "lon" in file '//lithfile(1:n))
    ierr = nf90_inq_varid( fileid, 'lat', varid )
    call isitok(ierr,'Error while getting "lat" variable identifiant in file '//lithfile(1:n))
    ierr = nf90_get_var( fileid, varid, ref_y_axis )
    call isitok(ierr,'Error while getting variable "lat" in file '//lithfile(1:n))

  ! File closing:
    ierr = nf90_close(fileid)

  end if

!=========================================================================!
! computation
!=========================================================================!

! climate interpolation

  if (CO2_level>0) then

    if ( CO2_level < co2climber(1) .or. CO2_level > co2climber(nclimber) ) then
      print *, 'ERROR: CO2 level cannot be outside the range of the climate fields'
      stop
    else

      call find_co2_interval( CO2_level, k1, k2, xi )

      do i=1,npixel
        if (areaclimber(i)>0) then
          Tclim(i) = (1-xi)*Tairclimber(i,k1) + xi*Tairclimber(i,k2)
          runclim(i) = (1-xi)*Runclimber(i,k1) + xi*Runclimber(i,k2)
            if (runclim(i).lt.0) runclim(i)=0
        end if
      end do

    end if

  end if


! dynsoil computation:

  xlevs = (/  ( dble(k)/dble(nDSlev) , k=nDSlev,1,-1 )  /)

  reg_thick     = fillvalue
  reg_x_surf    = fillvalue
  reg_tau_surf  = fillvalue
  reg_z_prof    = fillvalue
  reg_tau_prof  = fillvalue

  loc_minslope =  minval( slope, 1, slope>0 )

  loc_Esum = 0
  loc_Wsum = 0

  do i = 1,npixel
    if (areaclimber(i)>0) then

      ! put slope value on missing-points (Antarctica)
      if (slope(i)==fillvalue .or. slope(i)==0) slope(i)=loc_minslope

      do j = 1,nlitho

        reg_eros(j,i) = erosion( Tclim(i), runclim(i), slope(i) , j)
        reg_prod(j,i) = reg_prod_opt( Tclim(i), runclim(i), j )
        Kmain         = dissolution_constant( Tclim(i), runclim(i), j )

        if ( runclim(i) > 0 .and. reg_prod(j,i) > reg_eros(j,i) ) then
          ! EXP. HUMPED SPF
          ! n = 0
          ! reg_thick(j,i) = d1*log(reg_prod(j,i)/(knorm*reg_eros(j,i)))
          ! do while ( abs( ( reg_prod(j,i)*soil_prod_func(reg_thick(j,i) - reg_eros(j,i) )   /   reg_eros(j,i))    >    precis )
          !   reg_thick(j,i) = inversion(reg_thick(j,i),knorm*reg_eros(j,i)/reg_prod(j,i)
          !   n = n+1
          !   if (n==1000) then
          !     print *, 'inversion failed. ABORTED'
          !     stop
          !   end if
          ! end do
          ! EXPONENTIAL OR INVERSE SPF:
          reg_thick(j,i)    = eq_reg_thick( reg_prod(j,i), reg_eros(j,i) )
          reg_tau_surf(j,i) = reg_thick(j,i)/reg_eros(j,i)
          reg_x_surf(j,i)   = exp( -1*Kmain * ((reg_thick(j,i)/reg_eros(j,i))**(DS_sigma(j)+1)) / (DS_sigma(j)+1) )
          !
          reg_thick(j,i)    = reg_thick(j,i) / DS_scaling_factor()
          reg_tau_surf(j,i) = reg_tau_surf(j,i) / DS_scaling_factor()
        else
          reg_thick(j,i)  = 0
          reg_tau_surf(j,i)   = 0
          reg_x_surf(j,i) = 1
          reg_eros(j,i)   = 0 !I0 * (1-k1)/knorm
        end if

        reg_z_prof(1,j,i) = 0.
        reg_tau_prof(1,j,i) = 0.
        k = 1
        do while (k<nDSlev)
          if (xlevs(k+1) > reg_x_surf(j,i)+DS_epsl()) then
            k = k + 1
            reg_tau_prof(k,j,i) = (-1*(DS_sigma(j)+1)*log(xlevs(k))/Kmain)**(1/(DS_sigma(j)+1)) / DS_scaling_factor()
            reg_z_prof(k,j,i)   = reg_eros(j,i)*reg_tau_prof(k,j,i)
          else
            k = nDSlev
          end if
        end do

        loc_Esum     = loc_Esum     + litho_frac(j,i) * areaclimber(i) * reg_eros(j,i)
        loc_vol_Wsum = loc_vol_Wsum + litho_frac(j,i) * areaclimber(i) * reg_eros(j,i)*(1-reg_x_surf(j,i))
        loc_Wsum     = loc_Wsum     + litho_frac(j,i) * areaclimber(i) * reg_eros(j,i)*(1-reg_x_surf(j,i))*CaMg_rock(j)

      end do
    end if
  end do

  print *, 'erosion flux (Mt/y):' ,           1e-9   *  2500 *  loc_Esum
  print *, 'CO2 consumption flux (Tmol/y):' , 1e-12  *  loc_Wsum
  print *, 'weathering / critical erosion-limited weathering (-): ' , loc_vol_Wsum/loc_Esum



!=========================================================================!
! Write output
!=========================================================================!

! Output file creation:
  ierr = nf90_create(ofname, NF90_CLOBBER, fileid)
  call isitok(ierr,'Error while output file creation')

! Dimensions:
  ierr = nf90_def_dim( fileid, 'lon',     nlon,   Odimids(1) )
  call isitok(ierr,'Error while defining dimension "lon"')
  ierr = nf90_def_dim( fileid, 'lat',     nlat,   Odimids(2) )
  call isitok(ierr,'Error while defining dimension "lat"')
  ierr = nf90_def_dim( fileid, 'litho',   nlitho, Odimids(3) )
  call isitok(ierr,'Error while defining dimension "litho"')
  ierr = nf90_def_dim( fileid, 'xlevels', nDSlev,   Odimids(4) )
  call isitok(ierr,'Error while defining dimension "x"')

! Variables:
  ierr = nf90_def_var(fileid, 'lon',           NF90_FLOAT, Odimids(1), Ovarids(1))
  call isitok(ierr,'Error while defining variable "lon"')
  ierr = nf90_def_var(fileid, 'lat',           NF90_FLOAT, Odimids(2), Ovarids(2))
  call isitok(ierr,'Error while defining variable "lat"')
  ierr = nf90_def_var(fileid, 'litho',         NF90_INT,   Odimids(3), Ovarids(3))
  call isitok(ierr,'Error while defining variable "x"')
  ierr = nf90_def_var(fileid, 'xlevels',       NF90_FLOAT, Odimids(4), Ovarids(4))
  call isitok(ierr,'Error while defining variable "x"')
  ierr = nf90_def_var(fileid, 'reg_thickness', NF90_FLOAT, Odimids(1:3), Ovarids(5))
  call isitok(ierr,'Error while defining variable "reg_thickness"')
  ierr = nf90_def_var(fileid, 'tau_surf',      NF90_FLOAT, Odimids(1:3), Ovarids(6))
  call isitok(ierr,'Error while defining variable "tau_surf"')
  ierr = nf90_def_var(fileid, 'x_P_surf',      NF90_FLOAT, Odimids(1:3), Ovarids(7))
  call isitok(ierr,'Error while defining variable "x_cat_surf"')
  ierr = nf90_def_var(fileid, 'z',             NF90_FLOAT, Odimids,      Ovarids(8))
  call isitok(ierr,'Error while defining variable "z"')
  ierr = nf90_def_var(fileid, 'tau',           NF90_FLOAT, Odimids,      Ovarids(9))
  call isitok(ierr,'Error while defining variable "tau"')

! Variables attributes:
  ! Names:
  ierr = nf90_put_att(fileid, Ovarids(1) , 'name', 'longitude')
  call isitok(ierr,'Error while putting variable "lon" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'name', 'latitude')
  call isitok(ierr,'Error while putting variable "lat" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(3) , 'name', 'litho')
  call isitok(ierr,'Error while putting variable "litho" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(4) , 'name', 'x_levels')
  call isitok(ierr,'Error while putting variable "x" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(5) , 'name', 'initial regolith thickness')
  call isitok(ierr,'Error while putting variable "reg_thickness" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(6) , 'name', 'initial initial surface tau')
  call isitok(ierr,'Error while putting variable "tau_surf" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(7) , 'name', 'initial surface cation abundance')
  call isitok(ierr,'Error while putting variable "x_cat_surf" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(8) , 'name', 'initial z profile')
  call isitok(ierr,'Error while putting variable "z" attribute "name"')
  ierr = nf90_put_att(fileid, Ovarids(9) , 'name', 'initial tau profile')
  call isitok(ierr,'Error while putting variable "tau" attribute "name"')
  ! Axis:
  ierr = nf90_put_att(fileid, Ovarids(1) , 'axis',      'X')
  call isitok(ierr,'Error while putting variable "lon" attribute "axis" ')
  ierr = nf90_put_att(fileid, Ovarids(1) , 'nav_model', 'Default grid')
  ierr = nf90_put_att(fileid, Ovarids(1) , 'long_name', 'Longitude')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'axis',      'Y')
  call isitok(ierr,'Error while putting variable "lat" attribute "axis" ')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'nav_model', 'Default grid')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'long_name', 'Latitude')
  ierr = nf90_put_att(fileid, Ovarids(3) , 'long_name', 'Lithology')
  ierr = nf90_put_att(fileid, Ovarids(4) , 'axis',      'Z')
  call isitok(ierr,'Error while putting variable "x" attribute "axis" ')
  ierr = nf90_put_att(fileid, Ovarids(4) , 'positive',  'down')
  ! Units:
  ierr = nf90_put_att(fileid, Ovarids(1) , 'units', 'degrees_east')
  call isitok(ierr,'Error while putting variable "lon" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(2) , 'units', 'degrees_north')
  call isitok(ierr,'Error while putting variable "lat" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(3) , 'units', '-')
  call isitok(ierr,'Error while putting variable "litho" attribute "-" ')
  ierr = nf90_put_att(fileid, Ovarids(4) , 'units', '-')
  call isitok(ierr,'Error while putting variable "x" attribute "-" ')
  ierr = nf90_put_att(fileid, Ovarids(5) , 'units', 'm')
  call isitok(ierr,'Error while putting variable "reg_thickness" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(6) , 'units', 'y')
  call isitok(ierr,'Error while putting variable "tau_surf" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(7) , 'units', '-')
  call isitok(ierr,'Error while putting variable "x_cat-surf" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(8) , 'units', 'm')
  call isitok(ierr,'Error while putting variable "z" attribute "units" ')
  ierr = nf90_put_att(fileid, Ovarids(9) , 'units', 'y')
  call isitok(ierr,'Error while putting variable "tau" attribute "units" ')
  ! Missing-values:
  ierr = nf90_put_att(fileid, Ovarids(5) , '_FillValue', fillvalue)
  call isitok(ierr,'Error while putting variable "reg_thickness" attribute "Missing-value"')
  ierr = nf90_put_att(fileid, Ovarids(6) , '_FillValue', fillvalue)
  call isitok(ierr,'Error while putting variable "tau_surf" attribute "Missing-value"')
  ierr = nf90_put_att(fileid, Ovarids(7) , '_FillValue', fillvalue)
  call isitok(ierr,'Error while putting variable "x_cat_surf" attribute "Missing-value"')
  ierr = nf90_put_att(fileid, Ovarids(8) , '_FillValue', fillvalue)
  call isitok(ierr,'Error while putting variable "z" attribute "Missing-value"')
  ierr = nf90_put_att(fileid, Ovarids(9) , '_FillValue', fillvalue)
  call isitok(ierr,'Error while putting variable "tau" attribute "Missing-value"')

  ierr = nf90_enddef(fileid)
  call isitok(ierr,'Error while end of definition')

! Writting:
  ierr = nf90_put_var(fileid, Ovarids(1) , ref_x_axis )
  call isitok(ierr,'Error while putting variable "lon"')
  ierr = nf90_put_var(fileid, Ovarids(2) , ref_y_axis )
  call isitok(ierr,'Error wile putting variable "lat"')
  ierr = nf90_put_var(fileid, Ovarids(3) , (/(i,i=1,nlitho)/) )
  call isitok(ierr,'Error while putting variable "litho"')
  ierr = nf90_put_var(fileid, Ovarids(4) , xlevs )
  call isitok(ierr,'Error while putting variable "x"')

  do j = 1,nlat
    ierr = nf90_put_var( fileid , Ovarids(5) , transpose(real(reg_thick( : , 1+(j-1)*nlon : j*nlon ))) , &
                                                start=(/1,j,1/) , count=(/nlon,1,nlitho/)            )
    call isitok(ierr,'Error wile putting variable "reg_thickness"')
  end do
  do j = 1,nlat
    ierr = nf90_put_var( fileid , Ovarids(6) , transpose(real(reg_tau_surf( : , 1+(j-1)*nlon : j*nlon ))) , &
                                                start=(/1,j,1/) , count=(/nlon,1,nlitho/)            )
    call isitok(ierr,'Error wile putting variable "tau_surf"')
  end do
  do j = 1,nlat
    ierr = nf90_put_var( fileid , Ovarids(7) , transpose(real(reg_x_surf( : , 1+(j-1)*nlon : j*nlon ))) , &
                                                start=(/1,j,1/) , count=(/nlon,1,nlitho/)            )
    call isitok(ierr,'Error wile putting variable "reg_x_surf"')
  end do
  do j = 1,nlat
    do k = 1,nlitho
      ierr = nf90_put_var( fileid , Ovarids(8) , transpose(real(reg_z_prof( : , k , 1+(j-1)*nlon : j*nlon ))) , &
                                                  start=(/1,j,k,1/) , count=(/nlon,1,1,nDSlev/)            )
      call isitok(ierr,'Error wile putting variable "z"')
    end do
  end do
  do j = 1,nlat
    do k = 1,nlitho
      ierr = nf90_put_var( fileid , Ovarids(9) , transpose(real(reg_tau_prof( : , k , 1+(j-1)*nlon : j*nlon ))) , &
                                                  start=(/1,j,k,1/) , count=(/nlon,1,1,nDSlev/)            )
      call isitok(ierr,'Error wile putting variable "tau"')
    end do
  end do

  ! Output file closing
  ierr = nf90_close(fileid)
  call isitok(ierr,'Error while output file closing')




  ! REPORT:
  print *
  print *
  print *, 'created file '//ofname
  print *




end program
