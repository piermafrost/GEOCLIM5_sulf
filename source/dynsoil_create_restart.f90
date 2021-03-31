module dynsoil_create_restart_mod
implicit none

contains

subroutine dynsoil_create_restart( output_path, run_name, restart_name, varname, varunits, varmissvalname, varmissval,   &
                                   xlevs,longi,latit,reg_thick,reg_x_surf,reg_tau_surf,reg_z_prof,reg_tau_prof, t,tstart )
use netcdf_io_module, only: create, close_file, def_dim, def_var, put_att_text, put_att_real, enddef, put_var_int1D, &
                            put_var_real1D, put_var_real3D, put_var_real4D
use netcdf
use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev, scaling_factor

  include 'coupler.inc'
  integer, parameter:: npxl = nlon*nlat
  character(len=*), intent(in):: output_path, run_name, restart_name
  character(len=*), intent(in), dimension(:):: varname, varunits, varmissvalname
  double precision, intent(in), dimension(:):: varmissval
  double precision, intent(in), dimension(nDSlev):: xlevs
  double precision, intent(in), dimension(nlon)::   longi
  double precision, intent(in), dimension(nlat)::   latit
  double precision, intent(in), dimension(nlitho,npxl)::         reg_thick, reg_x_surf, reg_tau_surf
  double precision, intent(in), dimension(nDSlev,nlitho,npxl):: reg_z_prof, reg_tau_prof
  double precision, intent(in):: t, tstart
  double precision, dimension(nlitho,npxl)::         loc_reg_thick, loc_reg_tau_surf
  double precision, dimension(nDSlev,nlitho,npxl):: loc_reg_z_prof, loc_reg_tau_prof
  character(len=200):: fname
  integer:: dimid(4), varid(15)
  integer:: fid, i, j, k

  !===================================== OUTPUT VARIABLES LIST: ======================================!
  !  OUTPUT: X, Y, litho, xlevs, t, area, litho_frac, slope, temp, runoff, h_soil, x_surf, tau_surf,  !
  !  var #:  1  2  3      4      5  6     7           8      9     10      11      12      13         !
  !          z, tau, Reg_prod, Reg_eros, reg_P_diss, reg_P_eros, x_surf_eros, x_P_mean, reg_mean_age, !
  !          14 15   16        17        18          19          20           21        22            !
  !          Li_Friv,  Li_Fsp,   Li_driv                                                              !
  !          23        24        25                                                                   !
  !===================================================================================================!

  ! divide variables by scaling factor (backward transformation):
  where (reg_thick/=varmissval(11))
    loc_reg_thick    = reg_thick / scaling_factor
  else where
    loc_reg_thick    = varmissval(11)
  end where
  where (reg_tau_surf/=varmissval(13))
    loc_reg_tau_surf = reg_tau_surf / scaling_factor
  else where
    loc_reg_tau_surf = varmissval(13)
  end where
  ! for some obscure reasons, the where loop sometimes crashes in the following cases:
  do j = 1,npxl
    do i = 1,nlitho
      do k = 1,nDSlev
        if (reg_z_prof(k,i,j) /= varmissval(14)) then
          loc_reg_z_prof(k,i,j)   = reg_z_prof(k,i,j) / scaling_factor
        else
          loc_reg_z_prof(k,i,j)   = varmissval(14)
        end if
        if (reg_tau_prof(k,i,j) /= varmissval(15)) then
          loc_reg_tau_prof(k,i,j)   = reg_tau_prof(k,i,j) / scaling_factor
        else
          loc_reg_tau_prof(k,i,j)   = varmissval(15)
        end if
      end do
    end do
  end do
!  where (reg_z_prof/=varmissval(14))
!    loc_reg_z_prof   = reg_z_prof / scaling_factor
!  else where
!    loc_reg_z_prof   = varmissval(14)
!  end where
!  where (reg_tau_prof/=varmissval(15))
!    loc_reg_tau_prof = reg_tau_prof / scaling_factor
!  else where
!    loc_reg_tau_prof = varmissval(15)
!  end where


  ! output (restart) file name:
  fname = trim(output_path)//trim(restart_name)//trim(run_name)//'.nc'

  ! create output file
  call create( fname , fid )

  ! golbal attributes
  call put_att_text(fid, (/NF90_GLOBAL/), (/'title'/),      (/'DynSoil initial condition (restart)'/))
  call put_att_text(fid, (/NF90_GLOBAL/), (/'parent_run'/), (/trim(run_name)/))
  if (use_dynsoil_steady_state) then
    call put_att_text(fid, (/NF90_GLOBAL/), (/'DynSoil_mode'/),  (/'steady-state'/))
  else
    call put_att_text(fid, (/NF90_GLOBAL/), (/'DynSoil_mode'/),  (/'dynamic'/))
  end if
  call put_att_real(fid, (/NF90_GLOBAL/), (/'run_duration_yr'/), (/real(t)/))
  call put_att_real(fid, (/NF90_GLOBAL/), (/'run_time_yr'/),     (/real(t-tstart)/))

  ! define dimensions
  call def_dim( fid , varname(1:4), (/nlon,nlat,nlitho,nDSlev/) , dimid )

  ! define dimension variables
  call def_var( fid , varname(1:1), (/NF90_FLOAT/), dimid(1:1) , varid(1:1)  )
  call def_var( fid , varname(2:2), (/NF90_FLOAT/), dimid(2:2) , varid(2:2)  )
  call def_var( fid , varname(3:3), (/NF90_INT/),   dimid(3:3) , varid(3:3)  )
  call def_var( fid , varname(4:4), (/NF90_FLOAT/), dimid(4:4) , varid(4:4)  )

  ! define other variables
  call def_var( fid , varname(11:13), (/NF90_FLOAT/), dimid(1:3) , varid(11:13) ) ! h_soil, x_surf, tau_surf
  call def_var( fid , varname(14:15), (/NF90_FLOAT/), dimid      , varid(14:15) ) ! z_prof, tau_prof

  ! put dimension attributes:
  call put_att_text( fid, varid(1:1), (/'axis'/), (/'X'/)                    )
  call put_att_text( fid, varid(1:1), (/'nav_model'/), (/'Default grid'/)    )
  call put_att_text( fid, varid(2:2), (/'axis'/), (/'Y'/)                    )
  call put_att_text( fid, varid(2:2), (/'nav_model'/), (/'Default grid'/)    )
  call put_att_text( fid, varid(4:4), (/'axis'/), (/'Z'/)                    )
  call put_att_text( fid, varid(4:4), (/'positive'/), (/'down'/)             )

  ! put attributes
  call put_att_text( fid , (/varid(1:4),varid(11:15)/) ,    (/'name'/)        , (/varname(1:4),varname(11:15)/)     )
  call put_att_text( fid , (/varid(1:4),varid(11:15)/) ,    (/'units'/)       , (/varunits(1:4),varunits(11:15)/)   )
  call put_att_real( fid ,        varid(11:15)         , varmissvalname(11:15) , real(varmissval(11:15))            )

  ! end of definition
  call enddef( fid )

  ! put variables
  call put_var_real1D( fid, varid(1), real(longi)        )
  call put_var_real1D( fid, varid(2), real(latit)        )
  call put_var_int1D(  fid, varid(3), (/(j,j=1,nlitho)/) )
  call put_var_real1D( fid, varid(4), real(xlevs)        )
  call put_var_real3D( fid, varid(11), real(reshape(loc_reg_thick,   shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var_real3D( fid, varid(12), real(reshape(reg_x_surf,      shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var_real3D( fid, varid(13), real(reshape(loc_reg_tau_surf, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))) )
  call put_var_real4D( fid, varid(14), real(reshape(loc_reg_z_prof,   shape=(/nlon,nlat,nlitho,nDSlev/), order=(/4,3,1,2/))) )
  call put_var_real4D( fid, varid(15), real(reshape(loc_reg_tau_prof, shape=(/nlon,nlat,nlitho,nDSlev/), order=(/4,3,1,2/))) )

  ! close output file
  call close_file( fid )


end subroutine

end module
