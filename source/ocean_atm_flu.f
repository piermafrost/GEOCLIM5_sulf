    subroutine ocean_atm_flu(t)
!   ***************************
    use constante, only: dens
    implicit none
    include 'combine_foam.inc'

!   Thermodynamic constants:
    do j=1,nbasin - 1 ! pour ne pas calculer ds l atmos
        call eqcte(temp_box(j),salin(j),press_box(j),bco2(j) &
                   ,ak1(j),ak2(j),akb(j),akc(j))
!       carbonate speciation:
        call phbor(var_diss(1,j),var_diss(2,j),ch(j),salin(j),dens,ak1(j) &
                   ,ak2(j),akb(j))
        call chimie(var_diss(1,j),var_diss(2,j),ch(j),h2co3(j),hco3(j) &
                    ,co3(j),ak1(j),ak2(j))
        ph(j)=-dlog10(ch(j)*1.d-3)
        pco2_dissous(j)=h2co3(j)/bco2(j)
    end do
    return
    end
