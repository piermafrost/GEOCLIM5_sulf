    subroutine pco2_atm
!   *******************
    use constante, only: akk0, phias, phisa, PI_n_CO2_atm
    implicit none
    include 'combine_foam.inc'
    
    
    do j=1,nbasin
        fCO2atm_ocean(j)=akk0*(var_diss(7,nbasin)/PI_n_CO2_atm &
                         -pco2_dissous(j))*oce_surf(j) &
                         *indice_surface(j)
        fC13atm_ocean(j)=akk0*(phias*var_diss(7,nbasin)/PI_n_CO2_atm- &
                         (var_isot(1,j)-var_isot(4,nbasin)+phisa) &
                         *pco2_dissous(j))*oce_surf(j)*indice_surface(j)
    
        fC13ocean_atm(j)=akk0*((var_isot(4,nbasin)-var_isot(1,j)+phias) &
                         *var_diss(7,nbasin)/PI_n_CO2_atm-phisa*pco2_dissous(j)) &
                         *oce_surf(j)*indice_surface(j)
    
    enddo
    return
    end
