    subroutine diss_oxygen
!   **********************
    use constante, only: PI_n_O2_atm
    implicit none
    include 'combine_foam.inc'
    
!   full dynamic exchanges between surface and atmosphere
!   Calcul de Cste Henry pour oxygen
    do j=1,nbasin-1
        RT0=22.414d-03
        dum=temp_box(j)/100.
        xlogK=-58.3877+85.8079/dum+23.8439*dlog(dum) &
              +salin(j)*(-0.034892+dum*(-0.0019387*dum+0.015568))
        beta(j)=dexp(xlogK)/RT0
        beta(j)=beta(j)*0.21
    end do
    ! Set beta = 0 for atmospheric box to avoid numeric error while summing on all basins
    beta(nbasin) = 0

    

!   surface ocean at steady-state with the atmosphere:
    xnum=var_diss(6,nbasin)
    somme=0. 
    do j=1,nbasin
        somme = somme + beta(j)*box_vol(j)*indice_surface(j)
    enddo
    xden=PI_n_O2_atm + somme
    po2=xnum/xden   
    do j0=1,nsurface
        j = jbox_surface(j0)
        var_diss(6,j)=beta(j)*po2
    enddo
    return
    end
