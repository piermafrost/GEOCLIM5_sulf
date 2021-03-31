    subroutine sea_omega
!   ********************
    use constante, only: dens
    implicit none
    include 'combine_foam.inc'
    
    do j=1,nbasin-1
        omega(j)=var(4,j)*co3(j)/akc(j)
        omega_ara(j)=omega(j)/1.5
        call lyso(omega(j),temp_box(j),salin(j) &
                  ,dens,dplysc(j),dplysa(j))

    enddo
    return
    end
