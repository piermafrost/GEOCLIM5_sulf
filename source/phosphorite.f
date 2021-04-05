    subroutine phosphorite(t)
!   *************************
    use constante, only: akphos
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        fphos(i)=akphos*var(3,i)*clo*phosss*indice_deep(i)
        !!add O2 feedback
        !fphos(i)=akphos*var(3,i)*var(11,i)/2.161d-1*clo*phosss*indice_deep(i)
        !!          PO4 _^^^^^^^^ ^^^^^^^^^_ O2
    enddo
    return
    end
