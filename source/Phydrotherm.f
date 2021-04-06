    subroutine Phydrotherm(t)
!   *************************
!   scavenging of P on iron oxides (hydrothermal)
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        fhydP(i)=0.325d+10*var(3,i)/1.5d-3*clo*phosss* &
                 indice_deep(i)
        !! add O2 feedback
        !fhydP(i)=1.505d+10*var(3,i)/1.5d-3*var(11,i)*clo*phosss * indice_deep(i)
        !!              PO4 ^^^^^^^^     O2 ^^^^^^^^^
        !! stronger O2 feedback
        !fhydP(i)=6.965d+10*var(3,i)/1.5d-3*(var(11,i)**2) * clo*phosss * indice_deep(i)
        !!              PO4 ^^^^^^^^      O2 ^^^^^^^^^
    enddo
    return
    end
