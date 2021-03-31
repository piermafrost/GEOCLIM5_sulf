    subroutine salinity
!   *******************
    use constante, only: sal
    implicit none
    include 'combine_foam.inc'
    
    do j=1,nbasin
        salin(j)=sal
    enddo
    return
    end
