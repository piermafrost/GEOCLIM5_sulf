    subroutine phosphorite(t)
!   *************************
    use constante, only: akphos
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        fphos(i)=akphos*var(3,i)*clo*phosss*indice_deep(i)
    enddo
    return
    end
