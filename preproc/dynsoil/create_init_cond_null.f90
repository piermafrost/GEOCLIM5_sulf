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

program create_initial_condition_null

use netcdf
use isitiok_mod, only: isitok
implicit none


  include 'combine_foam.inc'
  character(len=*), parameter:: ofname='dynsoil_restart_nullreg.nc'
  real, parameter:: fillvalue = 9.96921e+36

  character(len=200):: area_fname
  double precision:: loc_area
  integer:: fileid, varid, ierr, Odimids(4), Ovarids(9)

!=========================================================================!

  ref_x_axis = (/ ( -180 + 360*(dble(i)-0.5)/nlon , i=1,nlon ) /)
  ref_y_axis = (/ ( -90  + 180*(dble(j)-0.5)/nlat , j=1,nlat ) /)
  xlevs = (/  ( k/dble(nDSlev) , k=nDSlev,1,-1 )  /)

!=========================================================================!


! INITIAL VALUE OF VARIABLES: 

  ! null regolith thickness:
  reg_thick = 0
  reg_x_surf = 1
  reg_tau_surf = 0
  reg_z_prof = fillvalue; reg_z_prof(1,:,:) = 0
  reg_tau_prof = reg_z_prof

!=========================================================================!


! MISSING-VALUE (area file argument received)

  n = command_argument_count()

  if (n>=1) then

    call get_command_argument( 1, area_fname )

    ! open area file
    open(unit=1,file=area_fname,status='old',action='read')

    do i=1,npixel

      ! load area
      read(unit=1,fmt=*) loc_area

      ! Put missing-value
      if (loc_area==0) then
        reg_thick(:,i)      = fillvalue
        reg_x_surf(:,i)     = fillvalue
        reg_tau_surf(:,i)   = fillvalue
        reg_z_prof(:,:,i)   = fillvalue
        reg_tau_prof(:,:,i) = fillvalue
      end if

    end do

    ! close area file
    close(unit=1)

  end if

!=========================================================================!

!!!!!!!!!!
! OUTPUT !
!!!!!!!!!!

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
  ierr = nf90_def_var(fileid, 'lon',       NF90_FLOAT, Odimids(1), Ovarids(1))
  call isitok(ierr,'Error while defining variable "lon"')
  ierr = nf90_def_var(fileid, 'lat',       NF90_FLOAT, Odimids(2), Ovarids(2))
  call isitok(ierr,'Error while defining variable "lat"')
  ierr = nf90_def_var(fileid, 'litho',     NF90_INT,   Odimids(3), Ovarids(3))
  call isitok(ierr,'Error while defining variable "x"')
  ierr = nf90_def_var(fileid, 'xlevels',   NF90_FLOAT, Odimids(4), Ovarids(4))
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




end program
