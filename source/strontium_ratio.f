    subroutine strontium_ratio(t)
!   *****************************
    use constante, only: rSrCar
    implicit none
    include 'combine_foam.inc'


    rsw=weighted_rsw
    do j0=1,nsurface
        j = jbox_surface(j0)
        rSrdep(j)=rSrCar*var(5,j)/80.d-3
    enddo
    return
    end
