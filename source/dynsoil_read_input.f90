!=================================================================================!
!------------------ SUBROUTINE FOR LOADING DYNSOIL INTPUT DATA -------------------!
!=================================================================================!


module dynsoil_read_input_mod
implicit none

contains

subroutine dynsoil_read_input( dynsoil_init_mode, IOcondID , &
              lon,lat,xlevs,reg_z_s,reg_x_s,reg_tau_s,reg_z,reg_tau,slope, missingpoints, slopemissval )

  use io_module, only: read_io_condition, load_netcdf_dimvar, load_netcdf, load_netcdf_3D, load_netcdf_4D, check_axis
  use dynsoil_create_init_condition, only: create_init_condition
  use utils, only: add_path, read_comment

  use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev, scaling_factor
  integer, parameter:: npxl = nlon*nlat
  character(len=*), intent(in):: dynsoil_init_mode
  integer, intent(in):: IOcondID
  double precision, intent(inout), dimension(nlon):: lon
  double precision, intent(inout), dimension(nlat):: lat
  double precision, intent(out), dimension(nDSlev):: xlevs
  double precision, intent(out), dimension(npxl):: slope
  double precision, intent(out), dimension(nlitho,npxl):: reg_z_s,reg_x_s,reg_tau_s
  double precision, intent(out), dimension(nDSlev,nlitho,npxl):: reg_z,reg_tau
  logical, intent(out), dimension(npxl):: missingpoints
  double precision, intent(out):: slopemissval
  double precision, dimension(nlon):: loc_lon
  double precision, dimension(nlat):: loc_lat
  !
  integer, parameter:: flg=200, vlg=50, ulg=30, mlg=20
  integer, parameter:: nvar=9
  character(len=flg):: restart_fname, slope_fname
  character(len=vlg):: var(nvar)
  character(len=ulg):: units(nvar)
  character(len=mlg):: missvalname(nvar)
  double precision  :: missval(nvar)
  !
  character(len=flg):: dummy
  integer:: i, nlev
  !
  integer, dimension(5):: ERROR_HANDLING_OPTION
  common /error/ ERROR_HANDLING_OPTION

  ! Default name and value of missing-values (if unspecified in text input files)
  character(len=*), parameter:: defmissvalname = '_FillValue'
  double precision, parameter:: defmissval = 9.96921e+36



  missvalname = defmissvalname
  missval = defmissval


  if (dynsoil_init_mode(1:8) == 'startup:') then

    ! Skip 9 uncommented lines (restart file info)
    do i = 1,9
      call read_comment(IocondID)
      read(unit=IOcondID, fmt=*)
    end do


  elseif(dynsoil_init_mode == 'restart') then

    print *
    print *
    print *, 'Read DynSoil conditions and load initialization files'
    print *, '-----------------------------------------------------'
    print *


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%   getting input file names:   %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    ! read input file name
    call read_comment(IOcondID)
    read(unit=IOcondID, fmt=*) restart_fname

    ! read variables attributes (#1 to #9):
    do i=1,3 ! x, y and z variable
      call read_comment(IOcondID)
      read(unit=IOcondID, fmt=*) dummy, var(i)
    end do
    do i=4,8 ! non-axis DynSoil restart variables
      call read_comment(IOcondID)
      read(unit=IOcondID, fmt=*) dummy, var(i), missvalname(i)
    end do

    ! Add GEOCLIM root path if needed
    call add_path(restart_fname)


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
    !%         loading data:        %!
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

    call load_netcdf_dimvar(restart_fname, varnameX=var(1),varnameY=var(2),varnameZ=var(3), varX=loc_lon,varY=loc_lat,varZ=xlevs)
    ! if reference axis do not already exist, create then
    if (all(lon==0d0) .and. all(lat==0d0)) then
      lon = loc_lon
      lat = loc_lat
    else ! otherwise, check axis consistency
      call check_axis('DynSoil initialization file', loc_lon, loc_lat, lon, lat, ERROR_HANDLING_OPTION(1))
    end if
    call load_netcdf_3D(restart_fname, var(4), reg_z_s,   units(4), missvalname(4), missval(4), expected_units='m')
    call load_netcdf_3D(restart_fname, var(5), reg_x_s,   units(5), missvalname(5), missval(5), expected_units='dimensionless')
    call load_netcdf_3D(restart_fname, var(6), reg_tau_s, units(6), missvalname(6), missval(6), expected_units='yr')
    call load_netcdf_4D(restart_fname, var(7), reg_z,     units(7), missvalname(7), missval(7), expected_units='m')
    call load_netcdf_4D(restart_fname, var(8), reg_tau,   units(8), missvalname(8), missval(8), expected_units='yr')


    ! Multiply by scaling factor
    ! --------------------------

    where (reg_z_s/=missval(4))
      reg_z_s   = reg_z_s * scaling_factor
    end where
    where (reg_tau_s/=missval(6))
      reg_tau_s = reg_tau_s * scaling_factor
    end where
    where (reg_z/=missval(7))
      reg_z = reg_z * scaling_factor
    end where
    where (reg_tau/=missval(8))
      reg_tau = reg_tau * scaling_factor
    end where


    ! Get missing points:
    ! -------------------

    where( reg_z_s(1,:)   == missval(4) &
      .or. reg_x_s(1,:)   == missval(5) &
      .or. reg_tau_s(1,:) == missval(6) &
      .or. reg_z(1,1,:)   == missval(7) &
      .or. reg_tau(1,1,:) == missval(8) )
       missingpoints = .true.
    elsewhere
       missingpoints = .false.
    end where


  else

    print *
    print *, 'ERROR: invalid DynSoil initialization mode "'//trim(dynsoil_init_mode)//'".'
    print *, 'Legal ones are "restart", "startup:null" and "startup:eq".'
    stop

  end if


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  !%             Slope            %!
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

  print *
  print *
  print *, 'Read Slope info and load data'
  print *, '-----------------------------'
  print *

  ! Read input file info
  call read_comment(IOcondID)
  read(unit=IOcondID, fmt=*) slope_fname, var(9), missvalname(9)

  ! Add GEOCLIM root path if needed
  call add_path(slope_fname)

  ! Load data
  if (all(lon==0d0) .and. all(lat==0d0)) then ! no reference axis
    call load_netcdf(slope_fname, var(9), slope, units(9), missvalname(9), missval(9), expected_units='dimensionless')
  else ! check consistency with reference axis
    call load_netcdf(slope_fname, var(9), slope, units(9), missvalname(9), missval(9), &
                     x_axis_ref=lon, y_axis_ref=lat, expected_units='dimensionless')
  end if
  slopemissval = missval(9)
  nlev = size(xlevs)


  ! Create init condition if "startup" mode
  ! ---------------------------------------

  if (dynsoil_init_mode(1:7)=='startup') then
    call create_init_condition(dynsoil_init_mode(9:), slopemissval, xlevs,reg_z_s,reg_x_s,reg_tau_s,reg_z,reg_tau, missingpoints)
    ! Notes:
    !   * remove "startup:" form 'dynsoil_init_mode' variable
    !   * scaling factor is already integrated in dynsoil physical parameters
    !     (for 'eq' initial conditions), and not needed for 'null' initial conditions
  end if



end subroutine

end module
