    subroutine oceanic_T(time)
!   +++++++++++++++++++++++++++++++
    use constante, only: PI_n_CO2_atm
    use co2_interpolation, only: find_co2_interval
    implicit none
    include 'combine_foam.inc'
    double precision:: xi

    ! Note: temp_box(nbasin) is the GMST (Global Mean Surface Temperature)

    if (nclimber==1) then ! fixed CO2 mode

        temp_box = Toceclimber(:,1)

    else ! CO2 interpolation

        p=var(12,nbasin)/PI_n_CO2_atm

        call find_co2_interval(co2climber, p, k1, k2, xi)
        temp_box  =  (1-xi) * Toceclimber(:,k1)  +  xi * Toceclimber(:,k2)

    end if

    end
