module dynsoil_steady_state
implicit none

contains


subroutine steady_state_surface_regolith(h,xs,taus, RP,E,FPdiss, temp,runoff,slope, dt, veget_factor,veget_eros_factor, &
                                         list_cont_pxl, ncontpxl)
use dynsoil_empirical_laws, only: erosion, reg_prod_opt, soil_prod_func, eq_reg_thick, dissolution_constant

use dynsoil_physical_parameters, only: nlon, nlat, nlitho, sigma, h0
integer, parameter:: npxl=nlon*nlat
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT-OUTPUT VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision, dimension(nlitho,npxl), intent(inout):: h, xs, taus                       !!
double precision, dimension(npxl),        intent(in):: temp, runoff, slope                  !!
double precision,                         intent(in):: dt                                   !!
!                                                                                           !!
double precision, dimension(nlitho,npxl), intent(out):: RP, E, FPdiss                       !!
! vegetation:                                                                               !!
double precision, dimension(npxl),        intent(in):: veget_factor, veget_eros_factor      !!
! continental points:                                                                       !!
integer, dimension(npxl), intent(in):: list_cont_pxl                                        !!
integer, intent(in):: ncontpxl                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision:: RPopt, Kmain, h2                                                !!
integer:: k, i, j0, j                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Continental loop
  ! ----------------
  do j0 = 1,ncontpxl
    j = list_cont_pxl(j0)

    ! test for non-null runoff points. If runoff is null, nothing is done
    ! because all fluxes (erosive, chemical...) are null.
    ! * * * * * * * * * * !
    if (runoff(j)>0) then !
    ! * * * * * * * * * * !

      ! Lithology loop
      ! --------------
      do i = 1,nlitho-1 ! Skip last lithology class => carbonates


        !====================================================!
        !========= climate-only dependent variables =========!
        !====================================================!

        ! optimal regolith production rate:
        RPopt  = veget_factor(j) * reg_prod_opt(temp(j), runoff(j), i)
        ! regolith erosion rate:
        E(i,j) = veget_eros_factor(j) * erosion(temp(j), runoff(j), slope(j), i)
        ! dissolution constant:
        Kmain  = veget_factor(j) * dissolution_constant(temp(j), runoff(j), i)


        !===============================================================================================================!
        !============== Equilibrium-corrected explicit Runge-Kutta 4 integration for regolith thickness (h) ============!
        !===============================================================================================================!
        !                                                                                                               !
        !                     Solve the equation  ' dh/dt = RPopt * soil_prod_func(h) - E '                             !
        !                                                                                                               !
        !===============================================================================================================!

        ! EXPONENTIAL SPF:

        ! QUASI STEADY-STATE: let regolith thickness evolve dynamically:
        ! analytical solution with constant RPopt and Erosion
        !h2   =   h(i,j)   + &
        !         h0 * log(  exp(-1*(E(i,j)/h0)*dt)  +  (RPopt/E(i,j)) * ( 1 - exp(-1*(E(i,j)/h0)*dt) ) * exp(-1*h(i,j)/h0) )
        !if (h2<0) h2 = 0
        !RP(i,j)   =   E(i,j)  +  ( h2 - h(i,j) ) / dt
        !h(i,j) = h2

        ! STEADY-STATE regolith thickness:
        h(i,j)    =   h0 * log(RPopt/E(i,j))
        RP(i,j)   =   E(i,j)
        if (h(i,j) < 0) h(i,j) = 0


        !===============================================================================================================!
        !==============                                    WEATHERING                                    ===============!
        !===============================================================================================================!
        !                                                                                                               !
        !           Solve the equations  ' dz/dt = RP + dissrate*dz/dx ' & ' dtau/dt = 1 + dissrate*dtau/dx '           !
        !           AT STEADY-STATE (West 2012 weathering model)                                                        !
        !                                                                                                               !
        !===============================================================================================================!


        taus(i,j) = h(i,j)/E(i,j)
        xs(i,j) = exp( -1*Kmain * ((h(i,j)/E(i,j))**(sigma(i)+1)) / (sigma(i)+1) )

        FPdiss(i,j) = RP(i,j)*(1-xs(i,j))


    end do
    ! end of lithology loop


    ! * * !
    else  !
    ! * * !

      ! null runoff: set all fluxes to 0
      h(1:nlitho-1, j)       = 0
      taus(1:nlitho-1, j)    = 0
      xs(1:nlitho-1, j)      = 1
      RP(1:nlitho-1, j)      = 0
      E(1:nlitho-1, j)       = 0
      FPdiss(1:nlitho-1, j)  = 0
    

    ! * * * !
    end if  !
    ! * * * !

  end do
  ! end of continental loop




end subroutine


subroutine steady_state_inner_regolith(x, h,xs,temp,runoff,E, z,tau,MPtot,ktop, veget_factor, &
                                       missingpoints)
use dynsoil_empirical_laws, only: dissolution_constant

use dynsoil_physical_parameters, only: nlon, nlat, nlitho, nDSlev,sigma, epsl 
integer, parameter:: npxl=nlon*nlat
!!
!!!!~!!!!!!!!!!!!!!!!!!!!!!!! INPUT-OUTPUT VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision, dimension(nDSlev),      intent(in):: x                        !!
double precision, dimension(nlitho,npxl), intent(in):: h, xs, E                 !!
double precision, dimension(npxl),        intent(in):: temp, runoff             !!
double precision, dimension(npxl),        intent(in):: veget_factor             !!
logical,          dimension(npxl),        intent(in):: missingpoints            !!
!                                                                               !!
double precision, dimension(nDSlev,nlitho,npxl), intent(out):: z, tau           !!
double precision, dimension(nlitho,npxl),        intent(out):: MPtot            !!
integer,          dimension(nlitho,npxl),        intent(out):: ktop             !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision:: mdx, Kmain                                                   !!
integer:: k, j, j0, i                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mdx = x(1)-x(2)

  ! Continental loop
  do j = 1,npxl
    if (.not. missingpoints(j)) then

      ! Lithology loop
      do i = 1,nlitho-1 ! Skip last lithology class => carbonates


        !====================================================!
        !========= climate-only dependent variables =========!
        !====================================================!

        ! dissolution constant:
        Kmain = veget_factor(j) * dissolution_constant(temp(j), runoff(j), i)


        !===============================================================================================================!
        !==============                            INNER REGOLITH COMPUTATION                            ===============!
        !===============================================================================================================!
        !                                                                                                               !
        !           Solve the equations  ' dz/dt = RP + dissrate*dz/dx ' & ' dtau/dt = 1 + dissrate*dtau/dx '           !
        !           AT STEADY-STATE (West 2012 weathering model)                                                        !
        !                                                                                                               !
        !===============================================================================================================!

        ktop(i,j)  = nDSlev
        MPtot(i,j) = 0
        z(1,i,j)   = 0
        tau(1,i,j) = 0
        k = 1
        do while ( k < ktop(i,j) )
          ! test if next point is outside the regolith:
          if ( x(k+1) < xs(i,j)+epsl ) then
            ktop(i,j) = k
            k = nDSlev
          else
            tau(k+1,i,j) = (-1*(sigma(i)+1) * log(x(k+1)) / Kmain)**(1/(sigma(i)+1))
            z(k+1,i,j)   = E(i,j)*tau(k+1,i,j)
            MPtot(i,j)   = MPtot(i,j)  +  mdx * (z(k+1,i,j) + z(k,i,j)) / 2
            k = k + 1
          end if
        end do
        MPtot(i,j) = MPtot(i,j)  +  (x(ktop(i,j)) - xs(i,j)) * (h(i,j) + z(ktop(i,j),i,j))/2   +   xs(i,j)*h(i,j)

      end do
    end if
  end do


end subroutine



end module
