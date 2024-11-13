module dynsoil
implicit none

contains



subroutine get_ktop_and_Pvol( x ,xsurf,h,z, missingpoints, ktop,reg_P_vol )
use dynsoil_physical_parameters, only: nlon, nlat, nDSlev, nlitho, epsl
integer, parameter:: npxl=nlon*nlat
double precision, intent(in), dimension(nDSlev):: x
double precision, intent(in), dimension(nlitho,npxl):: xsurf,h
double precision, intent(in), dimension(nDSlev,nlitho,npxl):: z
logical, intent(in), dimension(npxl):: missingpoints
integer, intent(out), dimension(nlitho,npxl):: ktop
double precision, intent(out), dimension(nlitho,npxl):: reg_P_vol
double precision:: mdx
integer:: i,j,k

mdx = x(1)-x(2)

! loop on all points:
do j = 1,npxl

  if (missingpoints(j)) then

    ktop(1:nlitho-1, j) = 0

  else

    ! loop on all lithologies expect carbonates:
    do i = 1,nlitho-1

      ktop(i,j) = nDSlev
      reg_P_vol(i,j) = 0
      k = 1
      do while ( k < nDSlev )
        ! test if next point is outside the regolith:
        if ( x(k+1) < xsurf(i,j)+epsl ) then
          ktop(i,j) = k
          k = nDSlev
        else
          reg_P_vol(i,j) = reg_P_vol(i,j)  +  mdx * (z(k+1,i,j) + z(k,i,j)) / 2
          k = k + 1
        end if
      end do
      reg_P_vol(i,j) = reg_P_vol(i,j)  +  (x(ktop(i,j)) - xsurf(i,j)) * (h(i,j) + z(ktop(i,j), i, j)) / 2  +  xsurf(i,j)*h(i,j)

    end do

  end if

end do


end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine regolith_time_integration( x, h,xs,taus,MPtot, ktop,z,tau, RP,E,FPdiss,FPeros,xs_eros, temp,runoff,slope, dt, &
                                      veget_factor, veget_eros_factor, list_cont_pxl, ncontpxl )

use dynsoil_empirical_laws, only: erosion, reg_prod_opt, soil_prod_func, eq_reg_thick, dissolution_constant
use dynsoil_physical_parameters, only: nlon, nlat, nDSlev, nlitho, h0, sigma, epsl

integer, parameter:: npxl = nlon*nlat
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT-OUTPUT VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision, dimension(nlitho,npxl),        intent(inout) :: h, xs, taus, MPtot              !!
integer,          dimension(nlitho,npxl),        intent(inout) :: ktop                            !!
double precision, dimension(nDSlev,nlitho,npxl), intent(inout) :: z, tau                          !!
!                                                                                                 !!
double precision, dimension(nDSlev),             intent(in) ::    x                               !!
double precision, dimension(npxl),               intent(in) ::    temp, runoff, slope             !!
double precision,                                intent(in)::     dt                              !!
!                                                                                                 !!
double precision, dimension(nlitho,npxl),        intent(out) ::   RP, E, FPdiss, FPeros, xs_eros  !!
! vegetation:                                                                                     !!
double precision, dimension(npxl),               intent(in) ::    veget_factor, veget_eros_factor !!
! continental points                                                                              !!
integer,          dimension(npxl),               intent(in) ::    list_cont_pxl                   !!
integer,                                         intent(in) ::    ncontpxl                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! INTERNAL VARIABLES: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
double precision:: heq, RPopt, Kmain, zsum, dissrate, Rdt, Ndt, O                  !!
double precision:: h1, h2, h3, RP1, RP2, RP3, RP4                                  !!
double precision:: mdx, xs2, taus2, z1(nDSlev), z2(nDSlev), tau2(nDSlev), MPtot2   !!
double precision:: z_noeros(nDSlev), zsum_noeros, MPtot_noeros                     !!
double precision:: xtop, ztop, tautop, taussigma_mean                              !!
double precision:: basalt_kin_boost                                                !!
integer:: ktopbis, ktop2                                                           !!
!                                                                                  !!
integer:: k, i, j, j0                                                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  mdx = x(1) - x(2)



  !######### Reminder: #########!
  !#  k -> vertical levels     #!
  !#  i -> lithology classes   #!
  !#  j -> continental points  #!
  !#############################!


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


        ! ! equilibrium estimated value:
        ! heq = eq_reg_thick( RPopt, E )
        !
        ! EXPONENTIAL HUMPED SPF:
        ! ! corrected RK4 scheme:
        ! if (h < heq) then
        !   ! correction so that h never goes beneath 0 and
        !   ! never exceed the (over)estimated value of h_equilibrium
        !   RP1 = RPopt*soil_prod_func(h)
        !   h1 = min(  heq  ,  max( 0. , h+(dt/2)*(RP1-E) )  )
        !   RP2 = RPopt*soil_prod_func(h1)
        !   h2 = min(  heq  ,  max( 0. , h+(dt/2)*(RP2-E) )  )
        !   RP3 = RPopt*soil_prod_func(h2)
        !   h3 = min(  heq  ,  max( 0. , h+dt*(RP3-E) )  )
        !   RP4 = RPopt*soil_prod_func(h3)
        !   RP = min(  (RP1+2*RP2+2*RP3+RP4)/6  ,  E + (heq-h)/dt  )
        ! else
        !   ! correction so that h never goes beneath 0
        !   RP1 = RPopt*soil_prod_func(h)
        !   h1 = max( 0. , h+(dt/2)*(RP1-E) )
        !   RP2 = RPopt*soil_prod_func(h1)
        !   h2 = max( 0. , h+(dt/2)*(RP2-E) )
        !   RP3 = RPopt*soil_prod_func(h2)
        !   h3 = max( 0. , h+dt*(RP3-E) )
        !   RP4 = RPopt*soil_prod_func(h3)
        !   RP = (RP1+2*RP2+2*RP3+RP4)/6
        ! end if
        ! E = min( E , h/dt + RP )
        ! h2 = h + dt*(RP-E)

        ! ! INVERSE SPF:
        ! ! euler implicit scheme:
        ! h2   =   ( h - E*dt + sqrt( (h-E*dt)**2 + 4*RPopt*h0*dt ) ) / 2
        ! RP = RPopt * soil_prod_func(h2)    

        ! EXPONENTIAL SPF:
        ! analytical solution with constant RPopt and Erosion
        h2        =   h(i,j)   + &
                      h0 * log(  exp(-1*(E(i,j)/h0)*dt)  +  (RPopt/E(i,j)) * ( 1 - exp(-1*(E(i,j)/h0)*dt) ) * exp(-1*h(i,j)/h0) )
        if (h2<0) h2 = 0
        RP(i,j)   =   E(i,j)  +  (h2 - h(i,j)) / dt




        !===============================================================================================================!
        !============== Euler t-implicit x-upstream integration for soil primary phases vertical profile ===============!
        !===============================================================================================================!
        !                                                                                                               !
        !           Solve the equations  ' dz/dt = RP + dissrate*dz/dx ' & ' dtau/dt = 1 + dissrate*dtau/dx '           !
        !           with dissrate = Kmain*x*tau^sigma (averaged between two x levels)                                   !
        !           + surface condition, where h(surf) is determined by reg prod/erosion, and x(surf) is computed       !
        !                                                                                                               !
        !===============================================================================================================!

        !###############################################  REMINDER  #########################################!
        !#                                                                                                  #!
        !#  nlev = number of vertical x levels                                                              #!
        !#  mdx = -dx = 1/nlev  { because x axis is decreasing from 1 to 0, and 0 is exclude from x axis }  #!
        !#  ktop = index of highest inner soil level                                                        #!
        !#  epsl = threshold for new point creation ( > 0 )                                                 #!
        !#                                                                                                  #!
        !#  *top refers to the highest * INSIDE THE REGOLTIH, and then a regular point of x axis            #!
        !#  whereas *s refers to * AT REGOLITH SURFACE, which position on x axis is 'xs' and height is 'h'  #!
        !#                                                                                                  #!
        !#  All computation are volumetric per unit of surface [m3/m2 = m]. Hence:                          #!
        !#  h      =  regolith thickness (m)                                                                #!
        !#             <=> volume of regolith per unit of surface                                           #!
        !#  MPtot  =  volume of primary phases per unit of surface (m)                                      #!
        !#             <=> equivalent thickness of primary phases                                           #!
        !#  RP     =  regolith production rate (m/y)                                                        #!
        !#             <=> volume of regolith created per year per unit of surface                          #!
        !#  E      =  regolith erosion rate (m/y)                                                           #!
        !#             <=> volume of regolith eroded per year per unit of surface                           #!
        !#  FPdiss =  volume of dissolved primary phases per year per unit of surface (m/y)                 #!
        !#             <=>* part of E affected by dissolution                                               #!
        !#  FPeros =  volume of eroded unaltered primary phases per year per unit of surface (m/y)          #!
        !#             <=>* part of E non affected by dissolution                                           #!
        !# * At steady-state only: RP = E = FPdiss + FPeros                                                 #!
        !#                                                                                                  #!
        !####################################################################################################!


        !*****************************************!
        ! INNER REGOLITH COMPUTATION (z and tau): !
        !*****************************************!

        !#### VERTICAL LOOP ####!
        !
        ! initialization
        z1(1) = 0
        z2(1) = 0
        tau2(1) = 0
        z_noeros(1) = 0
        zsum = 0
        zsum_noeros = 0
        ktopbis = 1
        !
        do k = 2,ktop(i,j)


          ! DISSOLUTION RATE, HEIGHTS (z) AND TIME (tau) OF x LEVELS:
          !----------------------------------------------------------
          !
          dissrate = averaged_dissolution_rate( Kmain, x(k-1), x(k), tau(k-1,i,j), tau(k,i,j), i )
          O = dissrate * dt/mdx
          tau2(k)  =  (  tau(k,i,j)  +          dt  +  O*tau2(k-1)  )  /  ( 1 + O )
          z1(k)    =  (  z(k,i,j)    +  RP(i,j)*dt  +  O*z1(k-1)    )  /  ( 1 + O )


          ! SUM OF z FOR MPtot COMPUTATION AND NO-EROSION COMPUTATION:
          !----------------------------------------------------------
          !
          if (h2 >= z1(k)) then ! still inside the regolith
          !
            ! NORMAL COMPUTATION
            ktopbis = k
            z_noeros(k) = z1(k)
            z2(k) = z1(k)
            zsum = zsum + z2(k)
            zsum_noeros = zsum
          !
          else           ! top of regolith reached
          !
            ! NO EROSION COMPUTATION:
            ! reactive dt:
            Rdt = dt *  (h(i,j) - z(k,i,j))  /  ( (h(i,j) - z(k,i,j))  -  (h2 - z1(k)) )
            ! non-reactive dt:
            Ndt = dt - Rdt
            !
            O = dissrate * Rdt/mdx
            z_noeros(k)   =   ( z(k,i,j)  +  RP(i,j)*Rdt +  O*z_noeros(k-1)   )  /  ( 1 + O )   +   RP(i,j)*Ndt
            z2(k) = h2
            zsum_noeros = zsum_noeros + z_noeros(k)
          !
          end if


        end do


        ! corrective term of the modified M sum:
        !---------------------------------------
        !
        MPtot2        =  mdx*zsum  +  z2(ktopbis) * (x(ktopbis) - mdx/2)
        MPtot_noeros  =  mdx*zsum_noeros  +  z_noeros(ktop(i,j)) * (x(ktop(i,j)) - mdx/2)



        !*********************************************!
        ! REGOLITH SURFACE COMPUTATION (xs and taus): !
        !*********************************************!


        if (h2 > 0) then


          ! dissolution rate at surface - x(surface) and tau(surface) :
          !------------------------------------------------------------
          !
          xtop   = x(ktopbis)
          ztop   = z2(ktopbis)
          tautop = tau2(ktopbis)
          !
          O = dt*E(i,j)/(h2-ztop)
          ! Euler implicit scheme for taus:
          taus2 = ( taus(i,j) + dt + O*tautop ) / ( 1 + O )
          ! Euler implicit scheme for xs:
          taussigma_mean = ( taus2**(sigma(i)+1) - tautop**(sigma(i)+1) ) / ((sigma(i)+1)*(taus2-tautop))
          xs2 = ( xs(i,j) + O*xtop ) / ( 1 + O + dt*Kmain*taussigma_mean)


          ! check for point addition:
          !--------------------------
          !
          k = ktopbis
          do while ( xs2 < x(k)-(1+epsl)*mdx )
            k = k+1
            if (k > ktop(i,j)) then
            ! point creation by linear interpolation :
              z2(k)   = ztop   + (h2-ztop)*(x(k)-xtop)/(xs2-xtop)
              tau2(k) = tautop + (taus2-tautop)*(x(k)-xtop)/(xs2-xtop)
            end if
          end do
          ktop2 = k


          ! Total mass (MPtot) and mass fluxes (FPdiss,FPeros):
          !-------------------------------------------------
          !
          if (ktopbis < ktop(i,j)) then
          ! one or several point(s) were removed
          !
            if (ktopbis < ktop2) then
              MPtot2 = MPtot2 + ((x(ktopbis+1)+x(ktopbis))/2)*(h2-z2(ktopbis))
            else
              MPtot2 = MPtot2 + ((xs2+x(ktop2))/2)*(h2-z2(ktop2))
            end if
            MPtot_noeros   =   MPtot_noeros   +   ((xs(i,j) + x(ktop(i,j)))/2) * ( (h(i,j) + RP(i,j)*dt)  -  z_noeros(ktop(i,j)) )
            FPeros(i,j)    =   (MPtot_noeros - MPtot2) / dt
            xs_eros(i,j)   =   FPeros(i,j) / E(i,j)
            FPdiss(i,j)    =   RP(i,j)  -  FPeros(i,j)  -  (MPtot2 - MPtot(i,j)) / dt
          !
          else
          ! no point removed
          !
            MPtot2        =  MPtot2  +  ((xs2 + x(ktop(i,j)))/2) * (h2 - ztop)
            xs_eros(i,j)  =  (xs(i,j) + xs2)/2
            FPeros(i,j)   =  xs_eros(i,j)*E(i,j)
            FPdiss(i,j)   =  RP(i,j)  -  FPeros(i,j)  -  (MPtot2 - MPtot(i,j)) / dt
          !
          end if


        else !(null regolith thickness)

          xs2          = 1
          taus2        = 0
          ktop2        = 1
          MPtot2       = 0
          xs_eros(i,j) = 1
          FPeros(i,j)  = E(i,j)
          FPdiss(i,j)  = 0

        end if




        !%%%%%%%%%%%%%%%%%%!
        ! VARIABLES UPDATE !
        !%%%%%%%%%%%%%%%%%%!
        !++++++++++++++++++!
        h(i,j)     = h2
        xs(i,j)    = xs2
        taus(i,j)  = taus2
        MPtot(i,j) = MPtot2
        ktop(i,j)  = ktop2
        z(:,i,j)   = z2
        tau(:,i,j) = tau2
        !++++++++++++++++++!


      end do
      ! end of litholoy loop


    ! * * !
    else  !
    ! * * !

      ! null runoff: set all fluxes to 0 for all lithology classes (expect carbonates)
      RP(1:nlitho-1, j)      = 0
      E(1:nlitho-1, j)       = 0
      FPdiss(1:nlitho-1, j)  = 0
      FPeros(1:nlitho-1, j)  = 0
      xs_eros(1:nlitho-1, j) = xs(1:nlitho-1, j)
    

    ! * * * !
    end if  !
    ! * * * !


  end do
  ! end of continental loop



end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


function averaged_dissolution_rate(K,x0,x1,tau0,tau1,klith)
! average the dissolution rate between two vertical x levels (x0 and x1)
! assuming a linear evolution of T between x0 and x1

  use dynsoil_physical_parameters, only: sigma
  double precision, intent(in):: K,x0, x1, tau0, tau1
  integer, intent(in):: klith
  double precision:: averaged_dissolution_rate

  ! poncutal dissolution rate to be averaged between x0 and x1:
  ! D(x) = K * x * tau(x)^sigma
  ! assuming that: tau(x) = tau0 + ((x-x0)/(x1-x0))*(tau1-tau0)
  averaged_dissolution_rate =   K * (                                      &
    x1*tau1**(sigma(klith)+1)  -  x0*tau0**(sigma(klith)+1)  -                           &
    (x1-x0) * (tau1**(sigma(klith)+2) - tau0**(sigma(klith)+2))/((sigma(klith)+2)*(tau1-tau0))  &
                                    )   /   ((sigma(klith)+1)*(tau1-tau0))

end function





end module
