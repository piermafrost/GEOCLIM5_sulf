subroutine printf(t,icompteur,y,ybio)
!*************************************
    use geoclim_write_output_mod, only: geoclim_write_output
    use geographic_write_output_mod, only: geographic_write_output
    use biodiv_write_output_mod,  only: biodiv_write_output
    use dynsoil_write_output_mod, only: dynsoil_write_output
    use dynsoil_offline_computation_mod, only: dynsoil_offline_computation
    use dynsoil_create_restart_mod, only: dynsoil_create_restart
    implicit none
    include 'combine_foam.inc'

!    compteur = compteur + 1
    icompteur = icompteur + 1
    if (float(icompteur).ge.xjump) then   !.or.isolver.eq.2) then
        icompteur=0
!        compteur=0

        !!!!! write geoclim outputs !!!!!
        call geoclim_write_output(t)

        !!!!! write biodiv outputs !!!!!
        if (coupling_ecogeo) then
            print *, 'ECOGEO module not available'
            stop
        end if

      
        !!!!! check for killing_signal !!!!!
        call read_deathnote(t)
        ! if killing signal is activated, tend <- t-10d6 and ageYprint <- t-10d6  =>  restart will be created, if not already done

    end if

    !!!!! write geographically-distributed outputs !!!!!
    icount_geogprint = icount_geogprint + 1
    if (icount_geogprint >= ijump_geogprint) then
        icount_geogprint = 0
        call geographic_write_output( GEOG_ofile_num, GEOG_ofile_name, GEOG_varout_name,  nlon, nlat, t, &
                                      Tclim, runclim, wth_allsil, wth_litho_wgh, wth_litho, fker, POC_export_rate, fp )
    end if


    !!!!! write dynsoil outputs (and create dynsoil offline variables): !!!!!
    if (coupling_dynsoil) then
        icount_DS_pri = icount_DS_pri + 1
        if ( icount_DS_pri >= ijump_DS_print ) then
            icount_DS_pri=0
            call dynsoil_offline_computation(xlevs, reg_z_prof,reg_tau_prof, reg_tau_surf,reg_x_surf,reg_thick,reg_eros, &
                                             Tclim,runclim,veget_factor, reg_ktop, reg_P_vol,                            &
                                             DS_missingpoints, DS_varout_missval(9),DS_varout_missval(10),               &
                                             reg_x_mean,reg_mean_age                                                     )
            call dynsoil_write_output(DS_ofile_num,DS_ofile_name,DS_varout_name, nlon,nlat, t, Tclim,runclim,  &
                                      reg_thick,reg_x_surf,reg_tau_surf, reg_z_prof,reg_tau_prof,              &
                                      reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros,                 &
                                      reg_x_mean,reg_mean_age, reg_Li_Friv,reg_Li_Fsp,reg_Li_driv              )
        end if
    end if


!!!! TEST FOR RESTART FILES: !!!!

    if (t.gt.ageYprint.and.fog.eq.0) then
        fog=1
        ! extensive variables (chemical amounts) => concentration
        do j = 1,11*nbasin
            write(10,*) y(j)/vol(j)
        end do
        ! PCO2 and isopotic variables
        do j = 1+11*nbasin, 19*nbasin
            write(10,*) y(j)
        end do
        ! rest of extensive variables => concentration
        do j = 1+19*nbasin,nvar
            write(10,*) y(j)/vol(j)
        end do
        ! ecological module restart
        if (coupling_ecogeo) then
            print *, 'ECOGEO module not available'
            stop
        end if
        ! DYNSOIL restart files:
        if (coupling_dynsoil)  then 
            call dynsoil_create_restart(output_path, run_name, DS_restart_name,                                                 &
                                        DS_varout_name, DS_varout_units, DS_varout_missvalname, DS_varout_missval,              &
                                        xlevs,ref_x_axis,ref_y_axis, reg_thick,reg_x_surf,reg_tau_surf,reg_z_prof,reg_tau_prof, &
                                        t, tsta                                                                                 )
        end if
    endif


!   do i=1,nbasin
!       bilan_water(i)=0.
!   enddo

!   do i=1,nbasin
!       do j=1,nbasin
!           bilan_water(i)=bilan_water(i)-F(i,j) &
!                          +F(j,i)
!       end do
!   enddo

!   do i=1,nbasin 
!       write(*,*)float(i),bilan_water(i)/(1.d+6*31.536d+6)
!   end do
!   stop


end subroutine
