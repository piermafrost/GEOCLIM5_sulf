    subroutine carb_dep(t)
!   **********************
    use constante, only: akcr, akdiss, hsurf, hthermo, hdeepepic, hdeep, rarag
    implicit none
    include 'combine_foam.inc'


    do j0=1,nappcont
        j = jbox_appcont(j0)
        if (omega_ara(j).ge.1.) then
            freef(j)=akcr*oce_surf(j)*1 &
             *(omega_ara(j)-1.)**1.7*clo                !Toar=18.33  !Ceno=24.16
            freefP(j)=(1./1000.)*freef(j)*clo*phosss*0. !Maas=16.66  !Berra=19.5
        else                                            !PT=22.0     !Carn=18.33
            freef(j)=0.                                 !Aptien=16.6 !Triasinf=13.33
            freefP(j)=0.                                !Rhetien=18.33
        endif
    end do



!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Open ocean Pelagic Carbonate dissolution (mol CaCO3/yr)

!   (area_** expressed in percentage of total seafloor)
!   no CCD correction !!

    do j=1,nbasin-1

        if (dplysc(j).lt.0.) dplysc(j)=0. !recycling efficiency cannot increase...
        if (dplysa(j).lt.0.) dplysa(j)=0. !...if lysocline rises above 0 m depth

        f_diss_c(j)=0.
        f_diss_a(j)=0.
        fdepc(j)=0.
        fdepa(j)=0.

    end do

    do j0=1,ndeep
        j = jbox_deep(j0)

        if (dplysc(j).lt.hdeep) then
             f_diss_c(j)=akdiss*(1.-dplysc(j) &
                         /hdeep)*(1.-rarag)
        else
            f_diss_c(j)=0.
        endif

        if (dplysa(j).lt.hdeep) then
            f_diss_a(j)=akdiss*(1.-dplysa(j) &
                        /hdeep)*rarag
        else
            f_diss_a(j)=0.
        endif

    end do



!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   thermocline reservoir: only recycling, no coastal interception :

    do j0=1,nthermo
        j = jbox_thermo(j0)

        if (dplysc(j).lt.hthermo) then
            f_diss_c(j)=akdiss*(1.-dplysc(j)/ &
                        hthermo)*(1.-rarag)
        else
            f_diss_c(j)=0.
        endif

        if (dplysa(j).lt.hthermo) then
            f_diss_a(j)=akdiss*(1.-dplysa(j)/ &
                        hthermo)*rarag
        else
            f_diss_a(j)=0.
        endif

    end do


!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   surface reservoir: only recycling, no coastal interception :

    do j0=1,nsurfnoappcont
        j = jbox_surfnoappcont(j0)

        if (dplysc(j).lt.hsurf) then
            f_diss_c(j)=akdiss*(1.-dplysc(j)/ &
                        hsurf)*(1.-rarag)
        else
            f_diss_c(j)=0.
        endif

        if (dplysa(j).lt.hsurf) then
            f_diss_a(j)=akdiss*(1.-dplysa(j)/ &
                        hsurf)*rarag
        else
            f_diss_a(j)=0.
        endif

    end do


!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   epicontinental surface:

    do j0=1,nappcont
        j = jbox_appcont(j0)

        if (dplysc(j).lt.hsurf) then
            f_diss_c(j)=akdiss*(1-dplysc(j)/hsurf) &
                        *(1.-rarag)
        else
            f_diss_c(j)=0.
        endif

        if (dplysa(j).lt.hsurf) then
            f_diss_a(j)=akdiss*(1-dplysa(j)/hsurf) &
                        *rarag
        else
            f_diss_a(j)=0.
        endif
    end do

!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   epicontinental deep:

    do j0=1,nnosurfbelappcont
        j = jbox_nosurfbelappcont(j0)

        if (dplysc(j).lt.hdeepepic) then
            f_diss_c(j)=akdiss*(1-dplysc(j)/hdeepepic) &
                        *(1.-rarag)
        else
            f_diss_c(j)=0.
        endif

        if (dplysa(j).lt.hdeepepic) then
            f_diss_a(j)=akdiss*(1-dplysa(j)/hdeepepic) &
                        *rarag
        else
            f_diss_a(j)=0.
        endif

    end do



    do j=1,nbasin

        fdissol_carb(j)=(f_diss_c(j)+f_diss_a(j))   !*shells ! fdownt
        fdissol_carbP(j)=(1./1000.)*fdissol_carb(j)*shells*0.
    enddo

    return
    end
