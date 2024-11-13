subroutine cont_weath(t)
!************************
use multidimensional_interpolation, only: climate_interpolation
use dynsoil, only: regolith_time_integration
use dynsoil_steady_state, only: steady_state_surface_regolith
use constante, only: PI_n_CO2_atm, Rgas, TSS_density, rsw_litho, CaMg_rock, P_rock, P2C_ker, P2C_carb, OC_rock, cp_cont, &
                     Sulf_rock, BASALT_LITHO_NUM
!use dynsoil_lithium, only: lithium
implicit none
double precision :: loc_granwth
include 'combine_foam.inc'
! WARNING: carbonate must be the last lithological class, whatever the chosen
! number of classes



! current atmospheric CO2 level (in PAL) and climate interpolation (from look-up table)
! -------------------------------------------------------------------------------------

p=var(12,nbasin)/PI_n_CO2_atm
call climate_interpolation(co2climber, clim_param_1, clim_param_2, clim_param_3, clim_param_4, clim_param_5,    &
                           p, cpvec, list_cont_pixel=list_cont_pixel, ncontpxl=ncontpxl,                        &
                           temp_array=Tairclimber, runf_array=Runclimber, interp_temp=Tclim, interp_runf=runclim)



! Total DISCHARGE of water to the ocean
!--------------------------------------

discharge=0.

do j0=1,ncontpxl
    j = list_cont_pixel(j0)
    discharge = discharge + 1e10*runclim(j)*areaclimber(j)
end do !                    ^^^^^
! convert area [1e6km2 -> m2] and runoff [cm/y -> m/y]






!========================================================!
!========================================================!
!==               WEATHERING COMPUTATION               ==!
!========================================================!
!========================================================!



call FACCO2(t,pco2,fco2,fe)





!=========================================================================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES DEPENDENT OF DYNSOIL USAGE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



if (.not. coupling_dynsoil) then
! ######################################
! # if not coupled with dynsoil module #
! ######################################


    ! initialising:
    fsilw=0.
    fbasw=0.
    fkerw=0.
    total_cont_POC_export=0.
    weighted_rsw=0.


    ! Riverine sediment discharge => parameteric power law of water discharge
    tss = 2.9d6 * discharge**0.5
     

    do j0=1,ncontpxl
        j = list_cont_pixel(j0)

        ! =====================
        ! Basaltic weathering
        ! =====================
        ! "basalt" lithology => number BASALT_LITHO_NUM
        facbas=(-42300./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        wth_litho(BASALT_LITHO_NUM,j) = (0.0483d+10)*dexp(facbas)*runclim(j) !0.044086200331d+10  !483685972.093689  !172421904
        wth_litho_wgh(BASALT_LITHO_NUM,j) = wth_litho(BASALT_LITHO_NUM,j) * litho_frac(BASALT_LITHO_NUM,j)
        fbas(j) = wth_litho_wgh(BASALT_LITHO_NUM,j) * areaclimber(j)
        fbasw=fbasw+fbas(j)*clo

        ! =====================
        ! Granitic weathering
        ! =====================
        ! skip basalts (BASALT_LITHO_NUM) and carbonates (last litho)
        facgra=(-48200./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        loc_granwth = (0.11379d+10)*dexp(facgra)*runclim(j) * fco2*fe !0.101357275d+10
        fsil(j) = 0
        do k = 1,BASALT_LITHO_NUM-1 ! skip basalts
            wth_litho(k,j) = loc_granwth
            wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
            fsil(j) = fsil(j) + wth_litho_wgh(k,j)
            weighted_rsw = weighted_rsw + rsw_litho(k)*wth_litho_wgh(k,j)*areaclimber(j)
        end do
        do k = BASALT_LITHO_NUM+1,nlitho-1 ! skip carbonates
            wth_litho(k,j) = loc_granwth
            wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
            fsil(j) = fsil(j) + wth_litho_wgh(k,j)
            weighted_rsw = weighted_rsw + rsw_litho(k)*wth_litho_wgh(k,j)*areaclimber(j)
        end do
        fsil(j) = fsil(j) * areaclimber(j)
        fsilw=fsilw+fsil(j)*clo

        ! conversion /1e6km2 => /m2 for ALL silicate weathering rates
        wth_litho(1:nlitho-1,j)     = 1e-12*wth_litho(1:nlitho-1,j)
        wth_litho_wgh(1:nlitho-1,j) = 1e-12*wth_litho_wgh(1:nlitho-1,j)

        wth_allsil(j) = sum(wth_litho_wgh(1:nlitho-1,j))

        ! =====================
        ! Kerogen weathering
        ! =====================
        ! Old formulation:
        !! facker=(-43000./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
        !fker(j)=4.4184d+8*runclim(j)*areaclimber(j)   !*var(11,nbasin)/38.d+18 !& !2.4 !1.432
        !                                              ! *dexp(facker)  !0.672d+10
        !fkerw=fkerw+fker(j)*clo
        fker(j) = 0.5 * sum(litho_frac(:,j)*OC_rock)*(tss/TSS_density/areatot)
        !         ^^^ : oxidation efficiency
        if (.not. lock_oxygen_cycle)  fkerw = fkerw + fker(j)*areaclimber(j)*clo
        fker(j) = 1e-12*fker(j) ! conversion /1e6km2 => /m2

    end do


    ! ======================
    ! Continental POC export
    ! ======================

    ! erosion export (t/km2/an)
    reg_eros_galy_unit=(tss/areatot)*1e-6*1e-3  !moving from kg/1e6km2/yr to t/km2/yr
    do j0=1,ncontpxl
        j = list_cont_pixel(j0)
        !scaling:
        POC_export_rate(j)=0.081*(reg_eros_galy_unit)**0.56  !in t/km2/yr
        POC_export_rate(j)=POC_export_rate(j)/12.         !mol/m2/yr
        total_cont_POC_export=total_cont_POC_export+POC_export_rate(j)*areaclimber(j)*1.e+12  !mol/yr
    end do


    ! Lithium riverine flux^M
    FrivLi=10.d+9

    ! lithology-averaged Strontium isotopic ratio
    weighted_rsw = weighted_rsw/fsilw


    ! ===============================================================================
    ! Pyrite weathering and distribution in carb dissol, sil dissol and H2SO4 release
    ! ===============================================================================

    ! molar flux of SO4^2- from pyrite oxidative weathering:
    fpyrw = 0.025*tss
    ! molar flux of Ca released to the ocean from carbonate dissolution by sulfuric acid
    fcarbsulfw  =  0.724*fpyrw
    ! molar flux of Ca released to the ocean from silicate dissolution by sulfuric acid
    fsilsulfw   =  0.276*fpyrw
    ! molar flux of H2SO4 from pyrite oxidation directly reaching the ocean
    fH2SO4sulfw =  0*fpyrw



else
! ######################################
! #   if coupled with dynsoil module   #
! ######################################


    ! increase time counter of dynsoil (=> asynchronous coupling) !
    icount_DS_int = icount_DS_int + 1
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


    if ( icount_DS_int >= ijump_DS_integration ) then ! if integration time has been reached

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! integrate dynsoil module !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! counter reinitialization:
        icount_DS_int = 0

        ! =============================
        ! DYNSOIL (silicate weathering)
        ! =============================
        !
        if (use_dynsoil_steady_state) then
            call steady_state_surface_regolith(reg_thick, reg_x_surf, reg_tau_surf, reg_prod, reg_eros, reg_P_diss, &
                                               Tclim, runclim, slope, DS_timestep, veget_factor, veget_eros_factor, &
                                               list_cont_pixel, ncontpxl                                            )
        else
            call regolith_time_integration(xlevs, reg_thick, reg_x_surf, reg_tau_surf, reg_P_vol, reg_ktop,                 &
                                           reg_z_prof, reg_tau_prof,                                                        &
                                           reg_prod,reg_eros,reg_P_diss,reg_P_eros, reg_x_surf_eros, Tclim, runclim, slope, &
                                           DS_timestep,                                                                     &
                                           veget_factor, veget_eros_factor, list_cont_pixel, ncontpxl                       )
        end if
        !call lithium(reg_P_diss, Tclim, runclim, reg_thick, reg_Li_Friv, reg_Li_Fsp, reg_Li_driv, list_cont_pixel, ncontpxl)


        ! Continental sum:
        ! ----------------

        ! initialising:
        fsilw=0.
        fbasw=0.
        FrivLi=0.
        dLiriv=0.
        tss=0.
        fkerw=0.
        total_cont_POC_export=0.
        weighted_rsw=0.
        fcarbsulfw=0
        fsilsulfw=0
        fH2SO4sulfw=0
        faddsulfw=0

        do j0=1,ncontpxl
            j = list_cont_pixel(j0)

            wth_allsil(j) = 0
            fsil(j) = 0
            do k = 1,nlitho-1 ! as k=nlitho is carbonate
                !
                wth_litho(k,j) = CaMg_rock(k) * reg_P_diss(k,j)
                wth_litho_wgh(k,j) = wth_litho(k,j) * litho_frac(k,j)
                wth_allsil(j) = wth_allsil(j) + wth_litho_wgh(k,j)
                fsil(j) = fsil(j)  +  wth_litho_wgh(k,j) * 1e12*areaclimber(j)
                !       WARNING: conversion: 1e6km2 -> m2: ^^^^
                weighted_rsw = weighted_rsw + rsw_litho(k)*wth_litho_wgh(k,j)*1e12*areaclimber(j)
                FrivLi = FrivLi + reg_Li_Friv(k,j)*litho_frac(k,j)*1e12*areaclimber(j)   !mol/yr
                dLiriv = dLiriv + reg_Li_driv(k,j)*(reg_Li_Friv(k,j)*litho_frac(k,j)*1e12*areaclimber(j))
            end do
            !
            ! Impose erosion of carbonates = erosion of sediments
            reg_eros(nlitho,j) = reg_eros(nlitho-1,j)
            !
            ! riverine sediment discharge
            tss = tss + TSS_density*sum(reg_eros(:,j)*litho_frac(:,j))*1e12*areaclimber(j)
            !
            ! Split silicates weathering in 2 fluxes: "basalt" weathering (mafic rocks) and "granite" weathering (rest of silicates)
            fbas(j) =  wth_litho_wgh(BASALT_LITHO_NUM,j) * 1e12*areaclimber(j)
            fsil(j) = fsil(j) - fbas(j)
            fsilw = fsilw + fsil(j)*clo
            fbasw = fbasw + fbas(j)*clo
            !
            ! same for strontium isotopic ratio
            weighted_rsw = weighted_rsw - rsw_litho(BASALT_LITHO_NUM)*wth_litho_wgh(BASALT_LITHO_NUM,j)*1e12*areaclimber(j)


            ! =====================
            ! Kerogen weathering
            ! =====================
            ! Old formulation:
            ! facker=(-43000./Rgas)*(1./(Tclim(j)+273.15)-1./288.15)
            ! fker(j)=4.4184d+8*runclim(j)*areaclimber(j)   !*var(11,nbasin)/38.d+18 !& !2.4 !1.432
                                                            ! *dexp(facker)  !0.672d+10
            fker(j) = 0.5 * sum(litho_frac(:,j)*OC_rock*reg_eros(:,j))
            !         ^^^ : oxidation efficiency
            if (.not. lock_oxygen_cycle)  fkerw = fkerw + fker(j)*1e12*areaclimber(j)*clo


            ! ======================
            ! Continental POC export
            ! ======================

            ! erosion export (t/km2/an)
            reg_eros_lithmean(j) = sum(litho_frac(:,j)*reg_eros(:,j)) / sum(litho_frac(:,j))
            !scaling:
            reg_eros_galy_unit=reg_eros_lithmean(j)*TSS_density*1.e6*1.e-3  !moving from m/yr to t/km2/yr
            POC_export_rate(j)=0.081*(reg_eros_galy_unit)**0.56  !in t/km2/yr
            POC_export_rate(j)=POC_export_rate(j)/12.         !mol/m2/yr
            total_cont_POC_export=total_cont_POC_export+POC_export_rate(j)*areaclimber(j)*1.e+12  !mol/yr


            ! ===============================================================================
            ! Pyrite weathering and distribution in carb dissol, sil dissol and H2SO4 release
            ! ===============================================================================

            ! molar flux of SO4^2- from pyrite oxidative weathering:
            fpyrw = sum(Sulf_rock(:) * reg_eros(:,j) * litho_frac(:,j))
            ! molar flux of Ca released to the ocean from carbonate dissolution by sulfuric acid
            fcarbsulfw  = fcarbsulfw    +  0.724*fpyrw * 1e12*areaclimber(j)
            ! molar flux of Ca released to the ocean from silicate dissolution by sulfuric acid
            fsilsulfw   = fsilsulfw     +  0.276*fpyrw * 1e12*areaclimber(j)
            ! molar flux of H2SO4 from pyrite oxidation directly reaching the ocean
            fH2SO4sulfw = fH2SO4sulfw   +  0*fpyrw  *1e12*areaclimber(j)

        end do

        ! Lithium riverine isotopic delta
        !dLiriv = dLiriv/FrivLi

        ! lithology-averaged Strontium isotopic ratio
        weighted_rsw = weighted_rsw/fsilw

    end if

end if




!=========================================================================!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VARIABLES INDEPENDENT OF DYNSOIL USAGE !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialising:
fcarbw=0.
fpw=0.
fpyrw=0.
ftrapw=0.

do j0=1,ncontpxl
    j = list_cont_pixel(j0)

    ! =====================
    ! Carbonate weathering
    ! =====================
    !fcarb(j)=0.18577793d+11*dsqrt(runclim(j))*areaclimber(j)*fco2*fe !old
    call carbo(p,Tclim(j),runclim(j),carb_weath_conc,j)
    ! Carbonate is assumed to be the last lithological class (k = nlitho)
    wth_litho(nlitho,j) = 3.589d10*runclim(j)*carb_weath_conc ! to fit 12.3 Tmol/yr of Carb wth
    !wth_litho(nlitho,j) = 16.*0.483483601042728d+10*runclim(j)*carb_weath_conc*0.978903 !full newton, old tuning
    wth_litho_wgh(nlitho,j) = litho_frac(nlitho,j) * wth_litho(nlitho,j)
    fcarb(j) = areaclimber(j) * wth_litho_wgh(nlitho,j)
    ! conversion because areaclimber is in 1e6km2 instead of m2:
    wth_litho(nlitho,j) = 1e-12 * wth_litho(nlitho,j)
    wth_litho_wgh(nlitho,j) = 1e-12 * wth_litho_wgh(nlitho,j)

    fcarbw=fcarbw+fcarb(j)*clo

    ! =====================
    ! Phosphorus weathering
    ! =====================
    fp(j) = sum((P_rock(1:nlitho-1)/CaMg_rock(1:nlitho-1)) * wth_litho_wgh(1:nlitho-1,j)) & ! silicate part
            +  P2C_carb * wth_litho_wgh(nlitho,j) & ! carbonate part
            +  P2C_ker * fker(j) ! kerogen part
    fpw = fpw + fp(j)*1d12*areaclimber(j)*clo*phosss

end do



! Substract Phosphorus that is in exported biospheric organic matter
fpw = fpw - total_cont_POC_export/cp_cont




! + + + + + + + + + !
! Lock sulfur cycle !
! + + + + + + + + + !
if (lock_sulfur_cycle) then
    ! => keep sil. sulf. wth (because it eventually consumes CO2), set H2SO4 release to 0 and
    !    use carb. sulf. wth as adjustment variable to balance sulfate-reduction
    fH2SO4sulfw = 0
    fcarbsulfw = sum( fSulfRed(1:nbasin-1) ) - fsilsulfw
end if


! + + + + + + + + + !
! Lock oxygen cycle !
! + + + + + + + + + !
if (lock_oxygen_cycle)  fkerw = sum( fO2_odc(1:nbasin-1) ) - (15./8.)*(fcarbsulfw + fsilsulfw + fH2SO4sulfw)






return
end

