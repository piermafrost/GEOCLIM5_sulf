    subroutine DOfA(O2,out)
!   ***********************
    implicit none
    include 'combine_foam.inc'

    out=4.7862d+01*O2*O2*O2*O2-3.6086d+1*O2*O2*O2 &
        +6.3448*O2*O2-2.2631*O2+1.0000

!   O.414 [] en moles / m3 de l'oxygene / degré d'anoxie

    if (O2.gt.0.414) then
        out=0.
    endif
!    ! no Oxygen P-feedback
!    out=0.495
    return
    end
