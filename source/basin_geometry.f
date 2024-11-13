    subroutine basin_geometry()
!    **************************
    use constante, only: ksed
    implicit none
    include 'combine_foam.inc'
    ! local variables:
    double precision:: area_sediepicontsurf, area_sediepicontnosurf, area_sedinoepicont

!   volumes in cubic meters (read in 1e6 km3):
    do k=1,nbasin
        read(33,*) box_vol(k)
        box_vol(k) = box_vol(k)*1.d+15 ! 1e6 km3 => m3
    end do

!   Define a volume for each variable (with special case for isotopic variables!)
!   COMBINE variable order: all dissolved, then all particulate, then all isotopic
    !
    i = 0
    do j=1,nvar_diss ! dissolved variables
        do k = 1,nbasin
            i = i+1
            vol(i) = box_vol(k)
        end do
    end do
    do j=1,nvar_part ! particulate variables
        do k = 1,nbasin
            i = i+1
            vol(i) = box_vol(k)
        end do
    end do
    do j=1,nvar_isot ! isotopic variables do not need volume => will be multiplied by 1
        do k = 1,nbasin
            i = i+1
            vol(i) = 1
        end do
    end do

    ndeep               = 0
    nnodeep             = 0
    nsedi               = 0
    nnosedi             = 0
    nthermo             = 0
    nnothermo           = 0
    nsurface            = 0
    nnosurface          = 0
    nepicont            = 0
    nnoepicont          = 0
    npolar              = 0
    nnopolar            = 0
    nappcont            = 0
    nnoappcont          = 0
    nsediepicontsurf    = 0
    nsediepicontnosurf  = 0
    nsedinoepicont      = 0
    nsurfappcont        = 0
    nsurfnoappcont      = 0
    nnosurfbelappcont   = 0
    nnosurfnobelappcont = 0

!   areas are read in 1e9 km2 and converted in m2
!   water fluxes are read in Sv (1e6 m3/s) and converted in m3/yr

    do i=1,nbasin

        read(34,*)oce_surf(i)
        read(35,*)indice_deep(i) ! =1 si oce profond sinon 0
        read(36,*)indice_sedi(i)
        read(37,*)surf_sedi(i)
        read(38,*)dummy, press_box(i)  !temperature is read in unit 32 "GCM_oceanic_temp.dat"
        read(39,*)app_cont(i)
        read(40,*)indice_thermo(i)
        read(41,*)indice_surface(i)
        read(42,*)(F(j,i),j=1,nbasin)
        read(44,*)fsink_inorg(i)
        read(45,*)fsink(i)
        read(46,*)indice_epicont(i)
        read(47,*)indice_polar(i)

        if(indice_deep(i)==1) then
            ndeep = ndeep + 1
            jbox_deep(ndeep) = i
        else
            nnodeep = nnodeep + 1
            jbox_nodeep(nnodeep) = i
        end if
        if(indice_sedi(i)==1) then
            nsedi = nsedi + 1
            jbox_sedi(nsedi) = i
            if (indice_epicont(i).eq.1) then
                if (indice_surface(i).eq.1) then
                    nsediepicontsurf = nsediepicontsurf + 1
                    jbox_sediepicontsurf(nsediepicontsurf) = i
                else
                    nsediepicontnosurf = nsediepicontnosurf + 1
                    jbox_sediepicontnosurf(nsediepicontnosurf) = i
                end if
            else
                nsedinoepicont = nsedinoepicont + 1
                jbox_sedinoepicont(nsedinoepicont) = i
            endif

        else
            nnosedi = nnosedi + 1
            jbox_nosedi(nnosedi) = i
        end if
        if(indice_thermo(i)==1) then
            nthermo = nthermo + 1
            jbox_thermo(nthermo) = i
        else
            nnothermo = nnothermo + 1
            jbox_nothermo(nnothermo) = i
        end if
        if(indice_surface(i)==1) then
            nsurface = nsurface + 1
            jbox_surface(nsurface) = i
            if(app_cont(i)==1) then
                nsurfappcont = nsurfappcont + 1
                jbox_surfappcont(nsurfappcont) = i
            else
                nsurfnoappcont = nsurfnoappcont + 1
                jbox_surfnoappcont(nsurfnoappcont) = i
            end if
        else
            nnosurface = nnosurface + 1
            jbox_nosurface(nnosurface) = i
            if(app_cont(i-1)==1) then
                nnosurfbelappcont = nnosurfbelappcont + 1
                jbox_nosurfbelappcont(nnosurfbelappcont) = i
            else
                nnosurfnobelappcont = nnosurfnobelappcont + 1
                jbox_nosurfnobelappcont(nnosurfnobelappcont) = i
            end if
        end if
        if(indice_epicont(i)==1) then
            nepicont = nepicont + 1
            jbox_epicont(nepicont) = i
        else
            nnoepicont = nnoepicont + 1
            jbox_noepicont(nnoepicont) = i
        end if
        if(indice_polar(i)==1) then
            npolar = npolar + 1
            jbox_polar(npolar) = i
        else
            nnopolar = nnopolar + 1
            jbox_nopolar(nnopolar) = i
        end if
        if(app_cont(i)==1) then
            nappcont = nappcont + 1
            jbox_appcont(nappcont) = i
        else
            nnoappcont = nnoappcont + 1
            jbox_noappcont(nnoappcont) = i
        end if

        oce_surf(i)=oce_surf(i)*1.d+15   ! 1e9 km2 => m2
        surf_sedi(i)=surf_sedi(i)*1.d+15 !
    enddo
    do i=1,nbasin
        do j=1,nbasin
            F(i,j)=F(i,j)*1.d+6*31.536d+6  ! Sv => m3/yr
        enddo
    enddo
    oce_surf_tot=oce_surf(nbasin)
! OBSOLETE
!    do n=1,nvar_real
!        read(43,*)indice_part(n)
!    enddo


    ! Sedimentation capacity of bottom bassins:
    ! -----------------------------------------

    ! Total sedimentation area
    area_sediepicontsurf = 0d0
    area_sediepicontnosurf = 0d0
    area_sedinoepicont = 0d0
    do j0 = 1,nsediepicontsurf
        j = jbox_sediepicontsurf(j0)
        area_sediepicontsurf = area_sediepicontsurf + surf_sedi(j) 
    end do
    do j0 = 1,nsediepicontnosurf
        j = jbox_sediepicontnosurf(j0)
        area_sediepicontnosurf = area_sediepicontnosurf + surf_sedi(j) 
    end do
    do j0 = 1,nsedinoepicont
        j = jbox_sedinoepicont(j0)
        area_sedinoepicont = area_sedinoepicont + surf_sedi(j) 
    end do

    ! fraction of total area of bassin type (-), and sedimentation capacity (in m3/yr)
    sedim_fract = 0d0
    sedim_capacity = 0d0
    do j0 = 1,nsediepicontsurf
        j = jbox_sediepicontsurf(j0)
        sedim_fract(j) = surf_sedi(j)/area_sediepicontsurf
        sedim_capacity(j) = sedim_fract(j) * ksed*area_sediepicontsurf**1.5
    end do
    do j0 = 1,nsediepicontnosurf
        j = jbox_sediepicontnosurf(j0)
        sedim_fract(j) = surf_sedi(j)/area_sediepicontnosurf
        sedim_capacity(j) = sedim_fract(j) * ksed*area_sediepicontnosurf**1.5
    end do
    do j0 = 1,nsedinoepicont
        j = jbox_sedinoepicont(j0)
        sedim_fract(j) = surf_sedi(j)/area_sedinoepicont
        sedim_capacity(j) = 1d21 ! infinity
    end do


    return
    end
