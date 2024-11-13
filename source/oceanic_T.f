    subroutine oceanic_T(time)
!   +++++++++++++++++++++++++++++++
    use constante, only: PI_n_CO2_atm
    use multidimensional_interpolation, only: climate_interpolation
    implicit none
    include 'combine_foam.inc'

    ! Note: temp_box(nbasin) is the GMST (Global Mean Surface Temperature)

    p = var_diss(7,nbasin)/PI_n_CO2_atm ! CO2 level in PAL

    call climate_interpolation(co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5, &
                               p, cpvec, boxtemp_array=Toceclimber, interp_boxtemp=temp_box                      )

    end
