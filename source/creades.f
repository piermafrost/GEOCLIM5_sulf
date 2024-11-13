    subroutine creades(x)
!   *********************

    use omp_lib
    use constante

    implicit none
    include 'combine_foam.inc'


    call oceanic_T(x)
    call biological_pump(x)
    call recycle(x)
    call salinity
    call ocean_atm_flu(x)
    call pco2_atm
    call diss_oxygen
    call sea_omega
    call carb_dep(x)
    call dc13_speciation
    call bio_frac
    call anoxic
    ! asynchronous coupling:
    if (icount_cont_weath==ijump_cont_weath) then
        icount_cont_weath = 0
        call cont_weath(x)
    else
        icount_cont_weath = icount_cont_weath + 1
    end if
    call org_dep(x)
    call degassing(x)
    call Phydrotherm(x)
    call phosphorite(x)
    call strontium_ratio(x)

    R_diss = 0.
    R_part = 0.
    R_isot = 0.

! WARNING: in order not to do the calculation in the atmospheric box, '-1' are
! added on some loop definition. It directly depends on indice_* files. Be
! carefull if you modified them. 




!<><><><><><><><><><><><><><>!
!<>  ===================   <>!
!<>  DISSOLVED VARIABLES   <>!
!<>  ===================   <>!
!<><><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! DIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=fCO2atm_ocean(j) - fbioC(j) &
                    - finorgC(j)+fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    +roxyd(j)*box_vol(j)*var_part(4,j) &
                    + fsink_inorg(j)*var_part(5,j)*(1.-closed)*box_vol(j) &
                    + fsink(j)*(box_vol(j))*var_part(4,j) &
                    *(1.-closed) &
                    +(fsink(j)*box_vol(j)*var_part(4,j)*indice_sedi(j)-fodc(j)) &
                    *clo &
                    + fmor(j) + fCO2_crust(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=fCO2atm_ocean(j) - fbioC(j) &
                    - finorgC(j)+fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    +roxyd(j)*box_vol(j)*var_part(4,j) &
                    + gotoshelf*fsink_inorg(j)*var_part(5,j) &
                    *(1.-closed)*box_vol(j) &
                    + gotoshelf*fsink(j)*(box_vol(j))*var_part(4,j) &
                    *(1.-closed) &
                    +(gotoshelf*fsink(j)*(box_vol(j))*var_part(4,j)-fodc(j)) &
                    *clo &
                    +2.*fsilw+2.*fbasw+2.*fcarbw-freef(j) &
                    +fkerw + fCO2_crust(j) + fcarbsulfw ! add Carbonate weathering by sulphuric acid

    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 ! ALK
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-2.*finorgC(j) &
                    +2.*fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    + 2.*fsink_inorg(j)*var_part(5,j) &
                    *(1.-closed)*(box_vol(j))-2.*fSO4_basin(j)-2.*fSO4_crust(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-2.*finorgC(j) &
                    +2.*fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    + 2.*gotoshelf*fsink_inorg(j)*var_part(5,j) &
                    *(1.-closed)*(box_vol(j)) &
                    +2.*fsilw+2.*fcarbw+2.*fbasw-2.*freef(j)-2.*fSO4_basin(j)-2.*fSO4_crust(j)-2.*fH2SO4sulfw
                                                                                            ! add H2SO4 from Pyrite weathering
    end do
    ! add Alkalinity produced by sulfate reduction in early diagenesis (2 alk for 1 S):
    do j0=1,nsedi
        j = jbox_sedi(j0)
        R_diss(i,j) = R_diss(i,j) + 2.*fSulfRed(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 ! PO4
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-fbioP(j)-finorgP(j)+fdissol_carb(j) &
                    *var_part(3,j)*box_vol(j) &
                    +roxyd(j)*var_part(2,j)*box_vol(j) &
                    +fsink(j)*(box_vol(j))*var_part(2,j)*(1.-closed) &
                    +(fsink(j)*(box_vol(j))*var_part(2,j)*indice_sedi(j) &
                    -fodp(j))*clo &
                    -fhydP(j)-fphos(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-fbioP(j)-finorgP(j)+fdissol_carb(j) &
                    *var_part(3,j)*box_vol(j) &
                    +roxyd(j)*var_part(2,j)*box_vol(j) &
                    +gotoshelf*fsink(j)*(box_vol(j))*var_part(2,j)*(1.-closed) &
                    +(gotoshelf*fsink(j)*(box_vol(j))*var_part(2,j)-fodp(j)) &
                    *clo+fpw
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! Ca
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-finorgC(j)+fdissol_carb(j) &
                    *var_part(5,j)*box_vol(j) &
                    + fsink_inorg(j)*var_part(5,j)*(1.-closed)*(box_vol(j))

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-finorgC(j)+fdissol_carb(j) &
                    *var_part(5,j)*box_vol(j) &
                    + gotoshelf*fsink_inorg(j)*var_part(5,j) &
                    *(1.-closed)*(box_vol(j)) &
                    +fsilw+fbasw+fcarbw-freef(j) &
                    +fcarbsulfw+fsilsulfw ! add carbonate and silicate weathering by sulphuric acid
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-finorgC(j)*rSrdep(j)+fdissol_carb(j) &
                    *var_part(1,j)*box_vol(j) &
                    + fsink_inorg(j)*(box_vol(j))*var_part(1,j)*(1.-closed)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_diss(i,j)=-finorgC(j)*rSrdep(j)+fdissol_carb(j) &
                    *var_part(1,j)*box_vol(j) &
                    + gotoshelf*fsink_inorg(j)*(box_vol(j))*var_part(1,j) &
                    *(1.-closed) &
                    +(fsilw+fbasw+fsilsulfw)*rSrsil+(fcarbw+fcarbsulfw)*rSrCar &! add carbonate and silicate weathering by sulphuric acid
                    -freef(j)*rSrdep(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=6 ! Oxy
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

!    bassins profonds

    do j0=1,nsurface
        j = jbox_surface(j0)
            R_diss(i,j)= 0.
    end do
    do j0=1,nnosurface-1!<- SKIP ATMOSPHERE BOX! 
        j = jbox_nosurface(j0)

        closed = clo*indice_sedi((j)) + 1*(1-indice_sedi((j)))

        R_diss(i,j)=-rO2C*roxyd(j)*box_vol(j) &
                    *var_part(4,j)-fsink(j)*box_vol(j) &
                    *var_part(4,j)*rO2C*(1.-closed) &
                    -(fsink(j)*box_vol(j)*indice_sedi(j) &
                    *var_part(4,j)*rO2C-fO2_odc(j))* &
                    clo
    end do

!    Atmosphere + ocean de surface a l equilibre

    do k=1,nbasin-1
        R_diss(i,nbasin)=R_diss(i,nbasin)+ rO2C*fbioC(k) &
                         -gotoshelf*fsink(k)*var_part(4,k)*box_vol(k)*rO2C &
                         *(1.-closed)*app_cont(k) &
                         -(gotoshelf*fsink(k)*var_part(4,k)*box_vol(k)*rO2C &
                         -fO2_odc(k))*clo*app_cont(k)-fkerw*app_cont(k)  &
                         +total_cont_POC_export*rO2C*app_cont(k) &
                         -(15./8.)*(fcarbsulfw+fsilsulfw+fH2SO4sulfw)*app_cont(k) ! add pyrite oxidation flux
    end do

    do k=1,nbasin-1
        do j0=1,nsurface
            j = jbox_surface(j0)
            R_diss(i,nbasin) = R_diss(i,nbasin) + F(k,j)*var_diss(i,k) - F(j,k)*var_diss(i,j)
        end do
    end do



!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=7 ! PCO2
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    if (nclimber==1) then ! fixed atm CO2 run
        R_diss(i,nbasin) = 0
    else
        do k=1,nbasin
            R_diss(i,nbasin)=R_diss(i,nbasin)-fCO2atm_ocean(k)
        end do
        R_diss(i,nbasin)=R_diss(i,nbasin)-2.*fsilw+fvol-fcarbw-2.*fbasw+ftrap+fanthros &
                         -total_cont_POC_export
    end if


!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc
    i=8 ! SO4^2-
!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc        

    ! Epicontinental box:
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)
        R_diss(i,j) = fcarbsulfw + fsilsulfw + fH2SO4sulfw ! pyrite oxidation S flux
    end do

    ! Sedimentary boxes
    do j0=1,nsedi
        j = jbox_sedi(j0)
        R_diss(i,j) = R_diss(i,j) - fSulfRed(j)
    end do




!<><><><><><><><><><><><><><><>!
!<>  =====================   <>!
!<>  PARTICULATE VARIABLES   <>!
!<>  =====================   <>!
!<><><><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! SrPIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R_part(i,j)= finorgC(j)*rSrdep(j)  &
                    - fdissol_carb(j)*var_part(1,j)*box_vol(j) &
                    + fsink_inorg(j-1)*(box_vol(j-1))*var_part(1,j-1) &
                    - fsink_inorg(j)*(box_vol(j))*var_part(1,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R_part(i,j)= finorgC(j)*rSrdep(j)  &
                    - fdissol_carb(j)*var_part(1,j)*box_vol(j) &
                    + (1.-gotoshelf)*fsink_inorg(j-1) &
                    *(box_vol(j-1))*var_part(1,j-1) &
                    - fsink_inorg(j)*(box_vol(j))*var_part(1,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R_part(i,j)= finorgC(j)*rSrdep(j)  &
                    - fdissol_carb(j)*var_part(1,j)*box_vol(j) &
                    - fsink_inorg(j)*(box_vol(j))*var_part(1,j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 !POP
!cccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R_part(i,j)= fbioP(j) - roxyd(j)*box_vol(j)*var_part(2,j) &
                    +fsink(j-1)*(box_vol(j-1))*var_part(2,j-1) &
                    -fsink(j)*(box_vol(j))*var_part(2,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R_part(i,j)= fbioP(j) - roxyd(j)*box_vol(j)*var_part(2,j) &
                    +(1.-gotoshelf)*fsink(j-1)*(box_vol(j-1))*var_part(2,j-1) &
                    -fsink(j)*(box_vol(j))*var_part(2,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R_part(i,j)= fbioP(j) - roxyd(j)*box_vol(j)*var_part(2,j) &
                    -fsink(j)*(box_vol(j))*var_part(2,j) + total_cont_POC_export*app_cont(j)/cp_cont
    end do




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 !PIP
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX! 
        j = jbox_nosurfnobelappcont(j0)
        R_part(i,j)= finorgP(j) - fdissol_carb(j)*var_part(3,j)*box_vol(j) &
                    +fsink_inorg(j-1)*(box_vol(j-1))*var_part(3,j-1) &
                    -fsink_inorg(j)*(box_vol(j))*var_part(3,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R_part(i,j)= finorgP(j) - fdissol_carb(j)*var_part(3,j)*box_vol(j) &
                    +(1.-gotoshelf)*fsink_inorg(j-1)*(box_vol(j-1))*var_part(3,j-1) &
                    -fsink_inorg(j)*(box_vol(j))*var_part(3,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R_part(i,j)= finorgP(j) - fdissol_carb(j)*var_part(3,j)*box_vol(j) &
                    -fsink_inorg(j)*(box_vol(j))*var_part(3,j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! POC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX! these are the deep reservoirs
        j = jbox_nosurfnobelappcont(j0)
        R_part(i,j)= fbioC(j) - roxyd(j)*box_vol(j)*var_part(4,j) &
                    + fsink(j-1)*(box_vol(j-1))*var_part(4,j-1) &
                    - fsink(j)*(box_vol(j))*var_part(4,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir (reservoir 7)
        j = jbox_nosurfbelappcont(j0)
        R_part(i,j)= fbioC(j) - roxyd(j)*box_vol(j)*var_part(4,j) &
                    + (1.-gotoshelf)*fsink(j-1)*(box_vol(j-1))*var_part(4,j-1) &
                    - fsink(j)*(box_vol(j))*var_part(4,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R_part(i,j)= fbioC(j) - roxyd(j)*box_vol(j)*var_part(4,j) &
                    - fsink(j)*(box_vol(j))*var_part(4,j)+total_cont_POC_export*app_cont(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! PIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R_part(i,j)= finorgC(j) - fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    + fsink_inorg(j-1)*var_part(5,j-1)*(box_vol(j-1)) &
                    - fsink_inorg(j)*var_part(5,j)*(box_vol(j))
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R_part(i,j)= finorgC(j) - fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    + (1.-gotoshelf)*fsink_inorg(j-1)*var_part(5,j-1)*(box_vol(j-1)) &
                    - fsink_inorg(j)*var_part(5,j)*(box_vol(j))
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R_part(i,j)= finorgC(j) - fdissol_carb(j)*var_part(5,j)*box_vol(j) &
                    - fsink_inorg(j)*var_part(5,j)*(box_vol(j))
    end do




!<><><><><><><><><><><><><>!
!<>  ==================  <>!
!<>  ISOTOPIC VARIABLES  <>!
!<>  ==================  <>!
!<><><><><><><><><><><><><>!


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! DIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_isot(i,j)=(fc13ocean_atm(j) &
                    -fbioC(j)*(dh2co3(j)-epsiC(j)-var_isot(i,j)) &
                    -finorgC(j)*(dco3(j)-var_isot(i,j)) &
                    +fdissol_carb(j)*var_part(5,j)*box_vol(j)*(var_isot(i+1,j)-var_isot(i,j)) &
                    +roxyd(j)*box_vol(j)*var_part(4,j)*(var_isot(i+2,j)-var_isot(i,j)) &
                    + fsink_inorg(j)*box_vol(j)*var_part(5,j)*(var_isot(i+1,j)-var_isot(i,j)) &
                    *(1.-closed) &
                    + fsink(j)*box_vol(j)*var_part(4,j)*(var_isot(i+2,j)-var_isot(i,j)) &
                    *(1.-closed)+ &
                    (fsink(j)*box_vol(j)*var_part(4,j)*indice_sedi(j)-fodc(j)) &
                    *clo &
                    *(var_isot(i+2,j)-var_isot(i,j))+fmor(j)*(dcmor-var_isot(i,j))) &
                    /(var_diss(1,j)*box_vol(j))
    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_isot(i,j)=(fc13ocean_atm(j) &
                    -fbioC(j)*(dh2co3(j)-epsiC(j)-var_isot(i,j)) &
                    -finorgC(j)*(dco3(j)-var_isot(i,j)) &
                    +fdissol_carb(j)*var_part(5,j)*box_vol(j)*(var_isot(i+1,j)-var_isot(i,j)) &
                    +roxyd(j)*box_vol(j)*var_part(4,j)*(var_isot(i+2,j)-var_isot(i,j)) &
                    + gotoshelf*fsink_inorg(j)*box_vol(j)*var_part(5,j) &
                    *(var_isot(i+1,j)-var_isot(i,j)) &
                    *(1.-closed) &
                    + gotoshelf* &
                    fsink(j)*box_vol(j)*var_part(4,j)*(var_isot(i+2,j)-var_isot(i,j)) &
                    *(1.-closed) &
                    +(gotoshelf*fsink(j)*box_vol(j)*var_part(4,j)-fodc(j)) &
                    *(var_isot(i+2,j)-var_isot(i,j))*clo &
                    +2.*fsilw*(var_isot(4,nbasin)-var_isot(i,j)) &
                    +2.*fbasw*(var_isot(4,nbasin)-var_isot(i,j)) &
                    +fcarbw*(var_isot(4,nbasin)-var_isot(i,j)) &
                    +(fcarbw+fcarbsulfw)*(dccarbw-var_isot(i,j)) & ! add Carbonate weathering by sulphuric acid
                    +fkerw*(dckerw-var_isot(i,j)) &
                    -freef(j)*(dco3(j)-var_isot(i,j))) &
                    /(var_diss(1,j)*box_vol(j))
    end do

    do j=1,nbasin-1
        do k=1,nbasin-1
            R_isot(i,j) = R_isot(i,j)  +  F(k,j)*var_diss(1,k)*(var_isot(i,k) - var_isot(i,j))/(var_diss(1,j)*box_vol(j))
        end do
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 ! PIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        if (var_part(5,j).gt.1.e-6) then
            R_isot(i,j)= (finorgC(j)*(dco3(j)-var_isot(i,j)) &
                        +fsink_inorg(j-1)*(box_vol(j-1))*var_part(5,j-1) &
                        *(var_isot(i,j-1)-var_isot(i,j)))/(var_part(5,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do
    do j0=1,nnosurfbelappcont
        j = jbox_nosurfbelappcont(j0)
        if (var_part(5,j).gt.1.e-6) then
            R_isot(i,j)= (finorgC(j)*(dco3(j)-var_isot(i,j)) &
                        +fsink_inorg(j-1)*(1.-gotoshelf) &
                        *(box_vol(j-1))*var_part(5,j-1) &
                        *(var_isot(i,j-1)-var_isot(i,j)))/(var_part(5,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        if (var_part(5,j).gt.1.e-6) then
            R_isot(i,j)= (finorgC(j)*(dco3(j)-var_isot(i,j))) &
                        /(var_part(5,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do


!    write(*,*)fkerw/1.d+12,(fodc(j)/1.d+12,j=1,10)
!    write(*,*)'*****************************'
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 ! POC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        if (var_part(4,j).gt.1.e-6) then
            R_isot(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var_isot(i,j)) &
                        +fsink(j-1)*box_vol(j-1)*var_part(4,j-1)*(var_isot(i,j-1)-var_isot(i,j)))/ &
                        (var_part(4,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do
    do j0=1,nnosurfbelappcont
        j = jbox_nosurfbelappcont(j0)
        if (var_part(4,j).gt.1.e-6) then
            R_isot(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var_isot(i,j)) &
                        +(1.-gotoshelf)*fsink(j-1)*box_vol(j-1)*var_part(4,j-1) &
                        *(var_isot(i,j-1)-var_isot(i,j)))/ &
                     (var_part(4,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        if (var_part(4,j).gt.1.e-6) then
            R_isot(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var_isot(i,j))  &
                        +total_cont_POC_export*app_cont(j)*(var_isot(4,nbasin)-epsiCont-var_isot(i,j)))/ &
                     (var_part(4,j)*box_vol(j))
        else
            R_isot(i,j)=0.
        endif
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=4 ! PCO2 dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do k=1,nbasin
        R_isot(i,nbasin)=R_isot(i,nbasin)-fC13atm_ocean(k) &
                                     /var_diss(7,nbasin)
    end do

    R_isot(i,nbasin)=R_isot(i,nbasin)+(fvol*(dcvol-var_isot(i,nbasin))+ &
                     ftrap*(dctrap-var_isot(i,nbasin))- &
                     total_cont_POC_export*(-epsiCont)) &
                     /var_diss(7,nbasin)




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=5 ! 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))
        R_isot(i,j)=(fdissol_carb(j)*var_part(1,j)*box_vol(j) &
                    *(var_isot(i+1,j)-var_isot(i,j))/(9.43+var_isot(i+1,j)) &
                    + fsink_inorg(j)*box_vol(j)*var_part(1,j)*(var_isot(i+1,j)-var_isot(i,j)) &
                    /(9.43+var_isot(i+1,j)) &
                    *(1.-closed)+fmor(j)*rSrmor*(rmor-var_isot(i,j)) &
                    /(9.43+rmor)) &
                    *(9.43+var_isot(i,j))/(var_diss(5,j)*box_vol(j))
    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R_isot(i,j)=(fdissol_carb(j)*var_part(1,j)*box_vol(j)*(var_isot(i+1,j)-var_isot(i,j)) &
                    /(9.43+var_isot(i+1,j)) &
                    +gotoshelf*fsink_inorg(j)*(box_vol(j))*var_part(1,j) &
                    *(var_isot(i+1,j)-var_isot(i,j))/(9.43+var_isot(i+1,j)) &
                    *(1.-closed) &
                    +fsilw*rSrSil*(rsw-var_isot(i,j))/(9.43+rsw) &
                    +fbasw*rSrSil*(rbas-var_isot(i,j))/(9.43+rbas) &
                    +fcarbw*rSrCar*(rcw-var_isot(i,j))/(9.43+rcw) &
                    +fcarbsulfw*rSrCar*(rcw-var_isot(i,j))/(9.43+rcw) & ! add carbonate weathering by sulphuric acid
                    +fsilsulfw*rSrSil*(rsw-var_isot(i,j))/(9.43+rsw)) & ! add silicate weathering by sulphuric acid
                    *(9.43+var_isot(i,j))/(var_diss(5,j)*box_vol(j))
    end do
    
    do j=1,nbasin-1
        do k=1,nbasin-1
            R_isot(i,j)=R_isot(i,j)+(F(k,j)*var_diss(5,k)* &
                        (var_isot(i,k)-var_isot(i,j))/(9.43+var_isot(i,k))) &
                        *(9.43+var_isot(i,j))/(var_diss(5,j)*box_vol(j))
        end do
    end do




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
        i=6 ! PIC 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
            j = jbox_nosurfnobelappcont(j0)
            if (var_part(1,j).gt.1.e-6) then
                    R_isot(i,j)= (finorgC(j)*rSrdep(j)*(var_isot(i-1,j)-var_isot(i,j)) &
                                /(9.43+var_isot(i-1,j)) &
                                +fsink_inorg(j-1)*(box_vol(j-1))*var_part(1,j-1)*(var_isot(i,j-1)-var_isot(i,j)) &
                                /(9.43+var_isot(i,j-1))) &
                                *(9.43+var_isot(i,j))/(var_part(1,j)*box_vol(j))
            else
                    R_isot(i,j)=0.
            endif
    end do
    do j0=1,nnosurfbelappcont
            j = jbox_nosurfbelappcont(j0)
            if (var_part(1,j).gt.1.e-6) then
                    R_isot(i,j)= (finorgC(j)*rSrdep(j)*(var_isot(i-1,j)-var_isot(i,j)) &
                                /(9.43+var_isot(i-1,j)) &
                                +(1.-gotoshelf)*fsink_inorg(j-1)*(box_vol(j-1)) &
                                *var_part(1,j-1)*(var_isot(i,j-1)-var_isot(i,j)) &
                                /(9.43+var_isot(i,j-1))) &
                                *(9.43+var_isot(i,j))/(var_part(1,j)*box_vol(j))
            else
                    R_isot(i,j)=0.
            endif
    end do
    do j0=1,nsurface
            j = jbox_surface(j0)
            if (var_part(1,j).gt.1.e-6) then
                    R_isot(i,j)= (finorgC(j)*rSrdep(j)*(var_isot(i-1,j)-var_isot(i,j)) &
                                /(9.43+var_isot(i-1,j))) &
                                *(9.43+var_isot(i,j))/(var_part(1,j)*box_vol(j))
            else
                    R_isot(i,j)=0.
            endif
    end do




    return
    end
