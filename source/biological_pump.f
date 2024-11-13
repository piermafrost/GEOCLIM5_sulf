    subroutine biological_pump(t)
!   *****************************
    use constante, only: cpred
    implicit none
    include 'combine_foam.inc'


    do j=1,nbasin

        ! Impose fraction of calcifying primary producter
        carb_ratio(j) = 0.3

        reff(j)=max(0d0, (pco2_dissous(j)-0.2)/(pco2_dissous(j)-0.1))

    end do

    ! Reduce productivity in polar basins
    do j0=1,npolar
        j = jbox_polar(j0)
        reff(j) = reff(j)/2.
        carb_ratio(j) = carb_ratio(j)/2
    end do

    ! Compute fluxes
    do j=1,nbasin

        fbioP(j)=0.
        do k=1,nbasin
            fbioP(j) = fbioP(j) + reff(j)*F(k,j)*var_diss(3,k)
        enddo

        xkill=1.  !mass extinction, obsolete

        fbioP(j)=(fbioP(j) + fpw*app_cont(j))*indice_surface(j)*xkill
        fbioC(j)=fbioP(j)*cpred

        if (omega(j).lt.1.0) then  !checking the saturation state of the ocean with respect to carbonates
            rC_Corg(j)=0.
        else
            rC_Corg(j)=carb_ratio(j)*(omega(j)-1)/(0.4+(omega(j)-1))
        endif

        ! Shelf-flag: calcifying organisms exist or not
        !   -> shelfal = 0.01 in non-epicontinental boxes if ishelfal==1 (i.e., calcifyers on shelves only)
        !      shelfal = 1    if ishelfal==0, and in epicontinental box regardless of ishelfal
        shelfal = 1 - 0.99*(1-indice_epicont(j))*ishelfal

        finorgC(j)=fbioC(j)*rC_Corg(j)*shells*shelfal
        finorgP(j)=0 ! => pas de P ds les coquilles ! !finorgC(j)/1000.

    end do

    return
    end
