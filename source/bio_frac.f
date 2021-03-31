    subroutine bio_frac
!   *******************
    implicit none
    include 'combine_foam.inc'

!       variable Ep (Kump and Arthur, 1999) (steawal3 file):
    do j=1,nbasin-1
        epsiC(j)=1e-3*(12.03*dlog10(1e3*h2co3(j))+1.19)
!        epsiC(j)=25.-((159.5*var(3,j)+38.39)/h2co3(j))
!       (Berner et al, 2000)
!        epsiC(j)=epsiC(j)/1000.+3.*((var(11,nbasin)/1.)-1)/1000.
!      constant fractionation
!        epsiC(j)=0.012
!        epsiC(j)=epsiC(j)/1.d+3
    enddo

    epsiCont=0.025

    return
    end
