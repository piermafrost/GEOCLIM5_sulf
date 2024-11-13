module dynsoil_create_restart_mod
implicit none

contains

subroutine dynsoil_create_restart(output_path, run_name, restart_name, outvar_info, DYNS_restart_dim, DYNS_restart_var, &
                                  xlevs, longi, latit, reg_thick, reg_x_surf, reg_tau_surf, reg_z_prof, reg_tau_prof, t, tstart)
use netcdf
use io_module, only: netcdf_output_var, DEFAULT_FILLVAL_NAME
use netcdf_io_module, only: create, close_file, def_dim, def_var, put_att, enddef, put_var
use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev, scaling_factor

  include 'coupler.inc'
  integer, parameter:: npxl = nlon*nlat
  character(len=*), intent(in):: output_path, run_name, restart_name
  type(netcdf_output_var), dimension(:), intent(in):: outvar_info
  character(len=*), intent(in):: DYNS_restart_dim(4), DYNS_restart_var(5)
  double precision, intent(in), dimension(nDSlev):: xlevs
  double precision, intent(in), dimension(nlon)::   longi
  double precision, intent(in), dimension(nlat)::   latit
  double precision, intent(in), dimension(nlitho,npxl)::         reg_thick, reg_x_surf, reg_tau_surf
  double precision, intent(in), dimension(nDSlev,nlitho,npxl):: reg_z_prof, reg_tau_prof
  double precision, intent(in):: t, tstart
  double precision, dimension(nlitho,npxl)::         loc_reg_thick, loc_reg_tau_surf
  double precision, dimension(nDSlev,nlitho,npxl):: loc_reg_z_prof, loc_reg_tau_prof
  character(len=200):: fname
  integer:: dimid(4), varid(10)
  integer:: fid, i, j, k

  ! divide variables by scaling factor (backward transformation):
  where (reg_thick/=outvar_info(6)%fillval)
    loc_reg_thick    = reg_thick / scaling_factor
  else where
    loc_reg_thick    = outvar_info(6)%fillval
  end where
  where (reg_tau_surf/=outvar_info(8)%fillval)
    loc_reg_tau_surf = reg_tau_surf / scaling_factor
  else where
    loc_reg_tau_surf = outvar_info(8)%fillval
  end where
  ! for some obscure reasons, the commented "where" loop sometimes crashes in the following cases:
  do j = 1,npxl
    do i = 1,nlitho
      do k = 1,nDSlev
        if (reg_z_prof(k,i,j) /= outvar_info(9)%fillval) then
          loc_reg_z_prof(k,i,j)   = reg_z_prof(k,i,j) / scaling_factor
        else
          loc_reg_z_prof(k,i,j)   = outvar_info(9)%fillval
        end if
        if (reg_tau_prof(k,i,j) /= outvar_info(10)%fillval) then
          loc_reg_tau_prof(k,i,j)   = reg_tau_prof(k,i,j) / scaling_factor
        else
          loc_reg_tau_prof(k,i,j)   = outvar_info(10)%fillval
        end if
      end do
    end do
  end do
!  where (reg_z_prof/=outvar_info(9)%fillval)
!    loc_reg_z_prof   = reg_z_prof / scaling_factor
!  else where
!    loc_reg_z_prof   = outvar_info(9)%fillval
!  end where
!  where (reg_tau_prof/=outvar_info(10)%fillval)
!    loc_reg_tau_prof = reg_tau_prof / scaling_factor
!  else where
!    loc_reg_tau_prof = outvar_info(10)%fillval
!  end where


  ! output (restart) file name:
  fname = trim(output_path)//trim(restart_name)//trim(run_name)//'.nc'

  ! create output file
  call create( fname , fid )

  ! golbal attributes
  call put_att(  fid, NF90_GLOBAL, 'title',           attribute_text='DynSoil initial condition (restart)')
  call put_att(  fid, NF90_GLOBAL, 'parent_run',      attribute_text=trim(run_name))
  if (use_dynsoil_steady_state) then
    call put_att(fid, NF90_GLOBAL, 'DynSoil_mode',    attribute_text='steady-state')
  else
    call put_att(fid, NF90_GLOBAL, 'DynSoil_mode',    attribute_text='dynamic')
  end if
  call put_att(  fid, NF90_GLOBAL, 'run_duration_yr', attribute_numeric=t)
  call put_att(  fid, NF90_GLOBAL, 'run_time_yr',     attribute_numeric=t-tstart)

  ! define dimensions
  call def_dim( fid , 'lon', nlon   , dimid(1) )
  call def_dim( fid , 'lat', nlat   , dimid(2) )
  call def_dim( fid , 'litho', nlitho , dimid(3) )
  call def_dim( fid , 'xlevs', nDSlev , dimid(4) )

  ! define dimension variables
  call def_var( fid , DYNS_restart_dim(1), 'float', dimid(1:1) , (/.true./) , varid(1)  )
  call def_var( fid , DYNS_restart_dim(2), 'float', dimid(2:2) , (/.true./) , varid(2)  )
  call def_var( fid , DYNS_restart_dim(3), 'int',   dimid(3:3) , (/.true./) , varid(3)  )
  call def_var( fid , DYNS_restart_dim(4), 'float', dimid(4:4) , (/.true./) , varid(4)  )

  ! define other variables
  call def_var( fid , DYNS_restart_var(1),  'float', dimid(1:3) , (/.true./) , varid(6)  ) ! h_soil
  call def_var( fid , DYNS_restart_var(2),  'float', dimid(1:3) , (/.true./) , varid(7)  ) ! x_surf
  call def_var( fid , DYNS_restart_var(3),  'float', dimid(1:3) , (/.true./) , varid(8)  ) ! tau_surf
  call def_var( fid , DYNS_restart_var(4),  'float', dimid      , (/.true./) , varid(9)  ) ! z_prof
  call def_var( fid , DYNS_restart_var(5),  'float', dimid      , (/.true./) , varid(10) ) ! tau_prof

  ! put dimension attributes:
  call put_att( fid, varid(1), 'axis',      attribute_text='X'               )
  call put_att( fid, varid(1), 'nav_model', attribute_text='Default grid'    )
  call put_att( fid, varid(2), 'axis',      attribute_text='Y'               )
  call put_att( fid, varid(2), 'nav_model', attribute_text='Default grid'    )
  call put_att( fid, varid(4), 'axis',      attribute_text='Z'               )
  call put_att( fid, varid(4), 'positive',  attribute_text='down'            )

  ! put attributes
  !do i = 1,4
  !  call put_att( fid , varid(i) , 'name'  , attribute_text=varname(i)    )
  !  call put_att( fid , varid(i) , 'units' , attribute_text=varunits(i)   )
  !end do
  do i = 6,10
    call put_att( fid , varid(i) , 'name'               , attribute_text=outvar_info(i)%vname                           )
    call put_att( fid , varid(i) , 'units'              , attribute_text=outvar_info(i)%units                           )
    call put_att( fid , varid(i) , DEFAULT_FILLVAL_NAME , attribute_numeric=outvar_info(i)%fillval, convert2type='real' )
  end do

  ! end of definition
  call enddef( fid )

  ! put variables
  call put_var( fid, varid(1),  var_real1D=real(longi)        )
  call put_var( fid, varid(2),  var_real1D=real(latit)        )
  call put_var( fid, varid(3),   var_int1D=(/(j,j=1,nlitho)/) )
  call put_var( fid, varid(4),  var_real1D=real(xlevs)        )
  call put_var( fid, varid(6),  var_real3D=real(reshape(loc_reg_thick,    shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var( fid, varid(7),  var_real3D=real(reshape(reg_x_surf,       shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var( fid, varid(8),  var_real3D=real(reshape(loc_reg_tau_surf, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var( fid, varid(9),  var_real4D=real(reshape(loc_reg_z_prof,   shape=(/nlon,nlat,nlitho,nDSlev/), order=(/4,3,1,2/))) )
  call put_var( fid, varid(10), var_real4D=real(reshape(loc_reg_tau_prof, shape=(/nlon,nlat,nlitho,nDSlev/), order=(/4,3,1,2/))) )

  ! close output file
  call close_file( fid )


end subroutine

end module
