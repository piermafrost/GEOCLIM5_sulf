    subroutine dc13_speciation
!   **************************
    implicit none
    include 'combine_foam.inc'
    
    do j=1,nbasin-1
        edb(j)=(24.12-9866/temp_box(j))*1.d-3 !H2CO3 -> HCO3 
        ebc(j)=((653.627/(temp_box(j)-233.45)**2)+0.22)*1.d-3 ! HCO3 -> CO3
        dco3(j)=(var(1,j)*var(13,j)-ebc(j)*hco3(j) &
                -ebc(j)*h2co3(j)-edb(j)*h2co3(j))/ &
                var(1,j)
        dhco3(j)=dco3(j)+ebc(j)
        dh2co3(j)=dhco3(j)+edb(j)
    end do
    return
    end
