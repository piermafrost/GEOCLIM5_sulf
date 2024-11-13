module dynsoil_offline_computation_mod
implicit none

contains


subroutine dynsoil_offline_computation(xlevs, z, tau, tausurf, xsurf, hsoil, erosion, temperature, runoff, veget_factor, &
                                       ktop, reg_P_vol, missingpoints, z_fillval, tau_fillval, x_mean, reg_mean_age)
  use dynsoil_steady_state, only: steady_state_inner_regolith
  use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev
  include 'coupler.inc'
  integer, parameter:: npxl=nlon*nlat
  double precision, intent(in), dimension(nDSlev):: xlevs
  double precision, intent(in), dimension(nlitho,npxl):: tausurf,xsurf,hsoil,erosion
  double precision, intent(inout), dimension(nDSlev,nlitho,npxl):: z, tau
  double precision, intent(in), dimension(npxl):: temperature, runoff, veget_factor
  integer, intent(inout), dimension(nlitho,npxl):: ktop
  double precision, intent(inout), dimension(nlitho,npxl):: reg_P_vol
  logical, intent(in), dimension(npxl):: missingpoints
  double precision, intent(in):: z_fillval, tau_fillval
  double precision, intent(out), dimension(nlitho,npxl):: x_mean,reg_mean_age
  integer:: j, k

  if (use_dynsoil_steady_state) then
    ! Compute inner regolith variables (offline, as not used for computing weathering rate)
    call steady_state_inner_regolith(xlevs, hsoil, xsurf, temperature, runoff, erosion, z, tau, &
                                     reg_P_vol, ktop, veget_factor, missingpoints               )
  end if
  call xP_average(reg_P_vol, hsoil, missingpoints, x_mean)
  call tau_average(xlevs, tau, tausurf, xsurf, ktop, reg_P_vol, missingpoints, reg_mean_age)

  ! Fillvalue for x points higher than regolith surface
  do j = 1,npxl
    if (.not. missingpoints(j)) then
      do k = 1,nlitho-1 ! Skip last lithology class => carbonates
        z(  ktop(k,j)+1:nDSlev, k, j) = z_fillval
        tau(ktop(k,j)+1:nDSlev, k, j) = tau_fillval
      end do
    end if
  end do

end subroutine



!----------------------------------------------------------------------------------------------------------------------------------!



subroutine tau_average(x, tau, tausurf, xsurf, ktop, reg_P_vol, missingpoints, tau_ave)
! Vertical primary-mass-averaged of tau in the regolith.
  use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev
  integer, parameter:: npxl = nlon*nlat
  double precision, intent(in), dimension(nDSlev):: x
  double precision, intent(in), dimension(nDSlev,nlitho,npxl):: tau
  double precision, intent(in), dimension(nlitho,npxl):: tausurf, xsurf, reg_P_vol
  integer, intent(in), dimension(nlitho,npxl):: ktop
  logical, intent(in), dimension(npxl):: missingpoints
  double precision, intent(out), dimension(nlitho,npxl):: tau_ave
  double precision:: mdx, tau_sum
  integer:: i,j,k

  mdx = x(1) - x(2)

  do i = 1,npxl
    if (.not. missingpoints(i)) then

      do j = 1,nlitho

        if (reg_P_vol(j,i)>0) then

          tau_sum = 0
          ! sum of regolith inner points:
          do k = 1,ktop(j,i)-1
            tau_sum    =    tau_sum   +   (                                           &
               x(k+1)*tau(k+1,j,i)**2  -  x(k)*tau(k,j,i)**2  +                       &
               mdx * (tau(k+1,j,i)**3-tau(k,j,i)**3) / (3*(tau(k+1,j,i)-tau(k,j,i)))  &
                                          )   /   (2*(tau(k+1,j,i)-tau(k,j,i)))
          end do
          ! term between ktop and surface:
          k = ktop(j,i)
          tau_sum    =    tau_sum   +   (                                                         &
             xsurf(j,i)*tausurf(j,i)**2  -  x(k)*tau(k,j,i)**2  +                                 &
             (x(k)-xsurf(j,i)) * (tausurf(j,i)**3-tau(k,j,i)**3) / (3*(tausurf(j,i)-tau(k,j,i)))  &
                                        )   /   (2*(tausurf(j,i)-tau(k,j,i)))

          tau_ave(j,i) = tau_sum / reg_P_vol(j,i)

        else

          tau_ave(j,i) = 0

        end if

      end do

    end if
  end do

end subroutine



!----------------------------------------------------------------------------------------------------------------------------------!



subroutine xP_average(reg_P_vol, h, missingpoints, xP_ave)
  use dynsoil_physical_parameters, only: nlon, nlat, nlitho
  integer, parameter:: npxl = nlon*nlat
  double precision, intent(in), dimension(nlitho,npxl):: reg_P_vol, h
  logical, intent(in), dimension(npxl):: missingpoints
  double precision, intent(out), dimension(nlitho,npxl):: xP_ave
  integer:: k
  !
  do k=1,nlitho
    where (.not. missingpoints)
      where (h(k,:) > 0)
        xP_ave(k,:) = reg_P_vol(k,:) / h(k,:)
      else where
        xP_ave(k,:) = 1
      end where
    end where
  end do
  !
end subroutine



!----------------------------------------------------------------------------------------------------------------------------------!



end module
