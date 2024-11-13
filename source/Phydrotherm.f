    subroutine Phydrotherm(t)
!   *************************
!   scavenging of P on iron oxides (hydrothermal)
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        fhydP(i)=0.325d+10*var(3,i)/1.5d-3*clo*phosss* &
                 indice_deep(i)
    enddo
    return
    end
