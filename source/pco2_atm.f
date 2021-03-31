    subroutine pco2_atm
!   *******************
    use constante, only: akk0, phias, phisa, PI_n_CO2_atm
    implicit none
    include 'combine_foam.inc'
    
    
    do j=1,nbasin
        fCO2atm_ocean(j)=akk0*(var(12,nbasin)/PI_n_CO2_atm &
                         -pco2_dissous(j))*oce_surf(j) &
                         *indice_surface(j)
        fC13atm_ocean(j)=akk0*(phias*var(12,nbasin)/PI_n_CO2_atm- &
                         (var(13,j)-var(16,nbasin)+phisa) &
                         *pco2_dissous(j))*oce_surf(j)*indice_surface(j)
    
        fC13ocean_atm(j)=akk0*((var(16,nbasin)-var(13,j)+phias) &
                         *var(12,nbasin)/PI_n_CO2_atm-phisa*pco2_dissous(j)) &
                         *oce_surf(j)*indice_surface(j)
    
    enddo
    return
    end
