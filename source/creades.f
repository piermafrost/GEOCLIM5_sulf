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
    call anoxic(x)
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

    do i=1,nvar_real
        do j=1,nbasin
            R(i,j)=0.
        end do
    end do


! WARNING: in order not to do the calculation in the atmospheric box, '-1' are
! added on some loop definition. It directly depends on indice_* files. Be
! carefull if you modified them. 

!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=1 ! DIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=fCO2atm_ocean(j) - fbioC(j) &
               - finorgC(j)+fdissol_carb(j)*var(10,j)*vol(j) &
               +roxyd(j)*vol(j)*var(9,j) &
               + fsink_inorg(j)*var(10,j)*(1.-closed)*vol(j) &
               + fsink(j)*(vol(j))*var(9,j) &
               *(1.-closed) &
               +(fsink(j)*vol(j)*var(9,j)*indice_sedi(j)-fodc(j)) &
               *clo &
               + fmor(j) + fCO2_crust(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=fCO2atm_ocean(j) - fbioC(j) &
               - finorgC(j)+fdissol_carb(j)*var(10,j)*vol(j) &
               +roxyd(j)*vol(j)*var(9,j) &
               + gotoshelf*fsink_inorg(j)*var(10,j) &
               *(1.-closed)*vol(j) &
               + gotoshelf*fsink(j)*(vol(j))*var(9,j) &
               *(1.-closed) &
               +(gotoshelf*fsink(j)*(vol(j))*var(9,j)-fodc(j)) &
               *clo &
               +2.*fsilw+2.*fbasw+2.*fcarbw-freef(j) &
               +fkerw + fCO2_crust(j)

    end do

!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=2 ! ALK
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-2.*finorgC(j) &
               +2.*fdissol_carb(j)*var(10,j)*vol(j) &
               + 2.*fsink_inorg(j)*var(10,j) &
               *(1.-closed)*(vol(j))-2.*fSO4_basin(j)-2.*fSO4_crust(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-2.*finorgC(j) &
               +2.*fdissol_carb(j)*var(10,j)*vol(j) &
               + 2.*gotoshelf*fsink_inorg(j)*var(10,j) &
               *(1.-closed)*(vol(j)) &
               +2.*fsilw+2.*fcarbw+2.*fbasw-2.*freef(j)-2.*fSO4_basin(j)-2.*fSO4_crust(j)-2.*fH2SO4sulfw
                                                                                       ! add H2SO4 from Pyrite weathering
    end do
    ! add Alkalinity produced by sulfate reduction in early diagenesis (2 alk for 1 S):
    do j0=1,nsedi
        j = jbox_sedi(j0)
        R(i,j) = R(i,j) + 2.*fSulfRed(j)
    end do



!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=3 ! PO4
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-fbioP(j)-finorgP(j)+fdissol_carb(j) &
               *var(8,j)*vol(j) &
               +roxyd(j)*var(7,j)*vol(j) &
               +fsink(j)*(vol(j))*var(7,j)*(1.-closed) &
               +(fsink(j)*(vol(j))*var(7,j)*indice_sedi(j) &
               -fodp(j))*clo &
               -fhydP(j)-fphos(j)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-fbioP(j)-finorgP(j)+fdissol_carb(j) &
               *var(8,j)*vol(j) &
               +roxyd(j)*var(7,j)*vol(j) &
               +gotoshelf*fsink(j)*(vol(j))*var(7,j)*(1.-closed) &
               +(gotoshelf*fsink(j)*(vol(j))*var(7,j)-fodp(j)) &
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

        R(i,j)=-finorgC(j)+fdissol_carb(j) &
               *var(10,j)*vol(j) &
               + fsink_inorg(j)*var(10,j)*(1.-closed)*(vol(j))

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-finorgC(j)+fdissol_carb(j) &
               *var(10,j)*vol(j) &
               + gotoshelf*fsink_inorg(j)*var(10,j) &
               *(1.-closed)*(vol(j)) &
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

        R(i,j)=-finorgC(j)*rSrdep(j)+fdissol_carb(j) &
               *var(6,j)*vol(j) &
               + fsink_inorg(j)*(vol(j))*var(6,j)*(1.-closed)

    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=-finorgC(j)*rSrdep(j)+fdissol_carb(j) &
               *var(6,j)*vol(j) &
               + gotoshelf*fsink_inorg(j)*(vol(j))*var(6,j) &
               *(1.-closed) &
               +(fsilw+fbasw+fsilsulfw)*rSrsil+(fcarbw+fcarbsulfw)*rSrCar &! add carbonate and silicate weathering by sulphuric acid
               -freef(j)*rSrdep(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=6 ! SrPIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R(i,j)= finorgC(j)*rSrdep(j)  &
               - fdissol_carb(j)*var(6,j)*vol(j) &
               + fsink_inorg(j-1)*(vol(j-1))*var(6,j-1) &
               - fsink_inorg(j)*(vol(j))*var(6,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R(i,j)= finorgC(j)*rSrdep(j)  &
               - fdissol_carb(j)*var(6,j)*vol(j) &
               + (1.-gotoshelf)*fsink_inorg(j-1) &
               *(vol(j-1))*var(6,j-1) &
               - fsink_inorg(j)*(vol(j))*var(6,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R(i,j)= finorgC(j)*rSrdep(j)  &
               - fdissol_carb(j)*var(6,j)*vol(j) &
               - fsink_inorg(j)*(vol(j))*var(6,j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=7 !POP
!cccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R(i,j)= fbioP(j) - roxyd(j)*vol(j)*var(7,j) &
               +fsink(j-1)*(vol(j-1))*var(7,j-1) &
               -fsink(j)*(vol(j))*var(7,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R(i,j)= fbioP(j) - roxyd(j)*vol(j)*var(7,j) &
               +(1.-gotoshelf)*fsink(j-1)*(vol(j-1))*var(7,j-1) &
               -fsink(j)*(vol(j))*var(7,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R(i,j)= fbioP(j) - roxyd(j)*vol(j)*var(7,j) &
               -fsink(j)*(vol(j))*var(7,j) + total_cont_POC_export*app_cont(j)/cp_cont
    end do




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=8 !PIP
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX! 
        j = jbox_nosurfnobelappcont(j0)
        R(i,j)= finorgP(j) - fdissol_carb(j)*var(8,j)*vol(j) &
               +fsink_inorg(j-1)*(vol(j-1))*var(8,j-1) &
               -fsink_inorg(j)*(vol(j))*var(8,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R(i,j)= finorgP(j) - fdissol_carb(j)*var(8,j)*vol(j) &
               +(1.-gotoshelf)*fsink_inorg(j-1)*(vol(j-1))*var(8,j-1) &
               -fsink_inorg(j)*(vol(j))*var(8,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R(i,j)= finorgP(j) - fdissol_carb(j)*var(8,j)*vol(j) &
               -fsink_inorg(j)*(vol(j))*var(8,j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=9 ! POC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX! these are the deep reservoirs
        j = jbox_nosurfnobelappcont(j0)
        R(i,j)= fbioC(j) - roxyd(j)*vol(j)*var(9,j) &
               + fsink(j-1)*(vol(j-1))*var(9,j-1) &
               - fsink(j)*(vol(j))*var(9,j)
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir (reservoir 7)
        j = jbox_nosurfbelappcont(j0)
        R(i,j)= fbioC(j) - roxyd(j)*vol(j)*var(9,j) &
               + (1.-gotoshelf)*fsink(j-1)*(vol(j-1))*var(9,j-1) &
               - fsink(j)*(vol(j))*var(9,j)
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R(i,j)= fbioC(j) - roxyd(j)*vol(j)*var(9,j) &
               - fsink(j)*(vol(j))*var(9,j)+total_cont_POC_export*app_cont(j)
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=10 ! PIC
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        R(i,j)= finorgC(j) - fdissol_carb(j)*var(10,j)*vol(j) &
               + fsink_inorg(j-1)*var(10,j-1)*(vol(j-1)) &
               - fsink_inorg(j)*var(10,j)*(vol(j))
    end do
    do j0=1,nnosurfbelappcont !we are below the epicontinental surface reservoir
        j = jbox_nosurfbelappcont(j0)
        R(i,j)= finorgC(j) - fdissol_carb(j)*var(10,j)*vol(j) &
               + (1.-gotoshelf)*fsink_inorg(j-1)*var(10,j-1)*(vol(j-1)) &
               - fsink_inorg(j)*var(10,j)*(vol(j))
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        R(i,j)= finorgC(j) - fdissol_carb(j)*var(10,j)*vol(j) &
               - fsink_inorg(j)*var(10,j)*(vol(j))
    end do

!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=11 ! Oxy
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

!    bassins profonds

    do j0=1,nsurface
        j = jbox_surface(j0)
            R(i,j)= 0.
    end do
    do j0=1,nnosurface-1!<- SKIP ATMOSPHERE BOX! 
        j = jbox_nosurface(j0)

        closed = clo*indice_sedi((j)) + 1*(1-indice_sedi((j)))

        R(i,j)=-rO2C*roxyd(j)*vol(j) &
               *var(9,j)-fsink(j)*vol(j) &
               *var(9,j)*rO2C*(1.-closed) &
               -(fsink(j)*vol(j)*indice_sedi(j) &
               *var(9,j)*rO2C-fO2_odc(j))* &
               clo
    end do

!    Atmosphere + ocean de surface a l equilibre

    do k=1,nbasin-1
        R(i,nbasin)=R(i,nbasin)+ rO2C*fbioC(k) &
                    -gotoshelf*fsink(k)*var(9,k)*vol(k)*rO2C &
                    *(1.-closed)*app_cont(k) &
                    -(gotoshelf*fsink(k)*var(9,k)*vol(k)*rO2C &
                    -fO2_odc(k))*clo*app_cont(k)-fkerw*app_cont(k)  &
                    +total_cont_POC_export*rO2C*app_cont(k) &
                    -(15./8.)*(fcarbsulfw+fsilsulfw+fH2SO4sulfw)*app_cont(k) ! add pyrite oxidation flux
    end do

    do k=1,nbasin-1
        do j0=1,nsurface
            j = jbox_surface(j0)
            R(i,nbasin)=R(i,nbasin)+F(k,j)*var(i,k) &
                        -F(j,k)*var(i,j)
        end do
    end do



!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=12 ! PCO2
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    if (nclimber==1) then ! fixed atm CO2 run
        R(i,nbasin) = 0
    else
        do k=1,nbasin
            R(i,nbasin)=R(i,nbasin)-fCO2atm_ocean(k)
        end do
        R(i,nbasin)=R(i,nbasin)-2.*fsilw+fvol-fcarbw-2.*fbasw+ftrap+fanthros &
                    -total_cont_POC_export + fcarbsulfw ! add Carbonate weathering by sulphuric acid
    end if


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=13 ! DIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=(fc13ocean_atm(j) &
               -fbioC(j)*(dh2co3(j)-epsiC(j)-var(i,j)) &
               -finorgC(j)*(dco3(j)-var(i,j)) &
               +fdissol_carb(j)*var(10,j)*vol(j)*(var(i+1,j)-var(i,j)) &
               +roxyd(j)*vol(j)*var(9,j)*(var(i+2,j)-var(i,j)) &
               + fsink_inorg(j)*vol(j)*var(10,j)*(var(i+1,j)-var(i,j)) &
               *(1.-closed) &
               + fsink(j)*vol(j)*var(9,j)*(var(i+2,j)-var(i,j)) &
               *(1.-closed)+ &
               (fsink(j)*vol(j)*var(9,j)*indice_sedi(j)-fodc(j)) &
               *clo &
               *(var(i+2,j)-var(i,j))+fmor(j)*(dcmor-var(i,j))) &
               /(var(1,j)*vol(j))
    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=(fc13ocean_atm(j) &
               -fbioC(j)*(dh2co3(j)-epsiC(j)-var(i,j)) &
               -finorgC(j)*(dco3(j)-var(i,j)) &
               +fdissol_carb(j)*var(10,j)*vol(j)*(var(i+1,j)-var(i,j)) &
               +roxyd(j)*vol(j)*var(9,j)*(var(i+2,j)-var(i,j)) &
               + gotoshelf*fsink_inorg(j)*vol(j)*var(10,j) &
               *(var(i+1,j)-var(i,j)) &
               *(1.-closed) &
               + gotoshelf* &
               fsink(j)*vol(j)*var(9,j)*(var(i+2,j)-var(i,j)) &
               *(1.-closed) &
               +(gotoshelf*fsink(j)*vol(j)*var(9,j)-fodc(j)) &
               *(var(i+2,j)-var(i,j))*clo &
               +2.*fsilw*(var(16,nbasin)-var(i,j)) &
               +2.*fbasw*(var(16,nbasin)-var(i,j)) &
               +fcarbw*(var(16,nbasin)-var(i,j)) &
               +fcarbw*(dccarbw-var(i,j)) &
               +fkerw*(dckerw-var(i,j)) &
               -freef(j)*(dco3(j)-var(i,j))) &
               /(var(1,j)*vol(j))
    end do

    do j=1,nbasin-1
        do k=1,nbasin-1
            R(i,j)=R(i,j)+F(k,j)*var(1,k)* &
                   (var(i,k)-var(i,j))/(var(1,j)*vol(j))
        end do
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=14 ! PIC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        if (var(10,j).gt.1.e-6) then
            R(i,j)= (finorgC(j)*(dco3(j)-var(i,j)) &
                   +fsink_inorg(j-1)*(vol(j-1))*var(10,j-1) &
                   *(var(i,j-1)-var(i,j)))/(var(10,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do
    do j0=1,nnosurfbelappcont
        j = jbox_nosurfbelappcont(j0)
        if (var(10,j).gt.1.e-6) then
            R(i,j)= (finorgC(j)*(dco3(j)-var(i,j)) &
                   +fsink_inorg(j-1)*(1.-gotoshelf) &
                   *(vol(j-1))*var(10,j-1) &
                   *(var(i,j-1)-var(i,j)))/(var(10,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        if (var(10,j).gt.1.e-6) then
            R(i,j)= (finorgC(j)*(dco3(j)-var(i,j))) &
                   /(var(10,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do


!    write(*,*)fkerw/1.d+12,(fodc(j)/1.d+12,j=1,10)
!    write(*,*)'*****************************'
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=15 ! POC dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_nosurfnobelappcont(j0)
        if (var(9,j).gt.1.e-6) then
            R(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var(i,j)) &
                   +fsink(j-1)*vol(j-1)*var(9,j-1)*(var(i,j-1)-var(i,j)))/ &
                   (var(9,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do
    do j0=1,nnosurfbelappcont
        j = jbox_nosurfbelappcont(j0)
        if (var(9,j).gt.1.e-6) then
            R(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var(i,j)) &
                   +(1.-gotoshelf)*fsink(j-1)*vol(j-1)*var(9,j-1) &
                   *(var(i,j-1)-var(i,j)))/ &
                (var(9,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do
    do j0=1,nsurface
        j = jbox_surface(j0)
        if (var(9,j).gt.1.e-6) then
            R(i,j)= (fbioC(j)*(dh2co3(j)-epsiC(j)-var(i,j))  &
                   +total_cont_POC_export*app_cont(j)*(var(16,nbasin)-epsiCont-var(i,j)))/ &
                (var(9,j)*vol(j))
        else
            R(i,j)=0.
        endif
    end do


!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=16 ! PCO2 dc13
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc

    do k=1,nbasin
        R(i,nbasin)=R(i,nbasin)-fC13atm_ocean(k) &
                                /var(12,nbasin)
    end do

    R(i,nbasin)=R(i,nbasin)+(fvol*(dcvol-var(i,nbasin))+ &
                ftrap*(dctrap-var(i,nbasin))- &
                total_cont_POC_export*(-epsiCont) &
                +fcarbsulfw*(dccarbw-var(i,nbasin))) & ! add Carbonate weathering by sulphuric acid
                /var(12,nbasin)




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    i=17 ! 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnoappcont
        j = jbox_noappcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=(fdissol_carb(j)*var(6,j)*vol(j) &
               *(var(i+1,j)-var(i,j))/(9.43+var(i+1,j)) &
               + fsink_inorg(j)*vol(j)*var(6,j)*(var(i+1,j)-var(i,j)) &
               /(9.43+var(i+1,j)) &
               *(1.-closed)+fmor(j)*rSrmor*(rmor-var(i,j)) &
               /(9.43+rmor)) &
               *(9.43+var(i,j))/(var(5,j)*vol(j))
    end do
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)

        closed = clo*indice_sedi(j) + 1*(1-indice_sedi(j))

        R(i,j)=(fdissol_carb(j)*var(6,j)*vol(j)*(var(i+1,j)-var(i,j)) &
               /(9.43+var(i+1,j)) &
               + gotoshelf*fsink_inorg(j)*(vol(j))*var(6,j) &
               *(var(i+1,j)-var(i,j))/(9.43+var(i+1,j)) &
               *(1.-closed) &
               +fsilw*rSrSil*(rsw-var(i,j))/(9.43+rsw) &
               +fbasw*rSrSil*(rbas-var(i,j))/(9.43+rbas) &
               +fcarbw*rSrCar*(rcw-var(i,j))/(9.43+rcw) &
               +fcarbsulfw*rSrCar*(rcw-var(i,j))/(9.43+rcw) & ! add carbonate weathering by sulphuric acid
               +fsilsulfw*rSrSil*(rsw-var(i,j))/(9.43+rsw)) & ! add silicate weathering by sulphuric acid
               *(9.43+var(i,j))/(var(5,j)*vol(j))
    end do

    do j=1,nbasin-1
        do k=1,nbasin-1
            R(i,j)=R(i,j)+(F(k,j)*var(5,k)* &
                   (var(i,k)-var(i,j))/(9.43+var(i,k))) &
                   *(9.43+var(i,j))/(var(5,j)*vol(j))
        end do
    end do




!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
        i=18 ! PIC 87Sr/86Sr
!cccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccc
    do j0=1,nnosurfnobelappcont-1!<- SKIP ATMOSPHERE BOX!
            j = jbox_nosurfnobelappcont(j0)
            if (var(6,j).gt.1.e-6) then
                    R(i,j)= (finorgC(j)*rSrdep(j)*(var(i-1,j)-var(i,j)) &
                           /(9.43+var(i-1,j)) &
                           +fsink_inorg(j-1)*(vol(j-1))*var(6,j-1)*(var(i,j-1)-var(i,j)) &
                           /(9.43+var(i,j-1))) &
                           *(9.43+var(i,j))/(var(6,j)*vol(j))
            else
                    R(i,j)=0.
            endif
    end do
    do j0=1,nnosurfbelappcont
            j = jbox_nosurfbelappcont(j0)
            if (var(6,j).gt.1.e-6) then
                    R(i,j)= (finorgC(j)*rSrdep(j)*(var(i-1,j)-var(i,j)) &
                           /(9.43+var(i-1,j)) &
                           +(1.-gotoshelf)*fsink_inorg(j-1)*(vol(j-1)) &
                           *var(6,j-1)*(var(i,j-1)-var(i,j)) &
                           /(9.43+var(i,j-1))) &
                           *(9.43+var(i,j))/(var(6,j)*vol(j))
            else
                    R(i,j)=0.
            endif
    end do
    do j0=1,nsurface
            j = jbox_surface(j0)
            if (var(6,j).gt.1.e-6) then
                    R(i,j)= (finorgC(j)*rSrdep(j)*(var(i-1,j)-var(i,j)) &
                           /(9.43+var(i-1,j))) &
                           *(9.43+var(i,j))/(var(6,j)*vol(j))
            else
                    R(i,j)=0.
            endif
    end do



!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc
    i=19 ! NOT ATTRIBUTED
!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc
    i=20 ! SO4^2-
!ccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccc        

    ! Epicontinental box:
    do j0=1,nappcont-1!<- SKIP ATMOSPHERE BOX!
        j = jbox_appcont(j0)
        R(i,j) = fcarbsulfw + fsilsulfw + fH2SO4sulfw ! pyrite oxidation S flux
    end do

    ! Sedimentary boxes
    do j0=1,nsedi
        j = jbox_sedi(j0)
        R(i,j) = R(i,j) - fSulfRed(j)
    end do




    return
    end
