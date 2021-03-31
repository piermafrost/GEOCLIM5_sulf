    subroutine biological_pump(t)
!   *****************************
    use constante, only: cpred
    implicit none
    include 'combine_foam.inc'
    

!   estimation of the ratio C:Corg in the organic productivity:
    a=0.3/0.2
    b=-1.5*a

    do j=1,nbasin
        if (omega(j).lt.1.0) then
            rC_Corg(j)=0.
        else
            rC_Corg(j)=0.3*(omega(j)-1.0)/(0.4+(omega(j)-1.))
        endif 
        reff(j)=max(0d0,(pco2_dissous(j)-0.2)/(pco2_dissous(j)-0.1))
    end do

    
    do j0=1,npolar
        j = jbox_polar(j0)
        reff(j)=reff(j)/2.
        if (omega(j).lt.1.0) then
            rC_Corg(j)=0.
        else
            rC_Corg(j)=0.15*(omega(j)-1.0)/(0.4+(omega(j)-1))
        endif 
    end do

    do j=1,nbasin
    
        fbioP(j)=0.
        do k=1,nbasin
            fbioP(j)= fbioP(j) + reff(j)*F(k,j)*var(3,k)
        enddo

        xkill=1.  !mass extinction, obsolete

        fbioP(j)=(fbioP(j) + fpw*app_cont(j))*indice_surface(j)*xkill


!       ***********************************************************
!       coupling GEOCLIM biodiv
        if (coupling_ecogeo) then
            print *, 'ECOGEO module not available'
            stop
        endif
!       ***********************************************************
            
        fbioC(j)=fbioP(j)*cpred
        
!       shelfal = 1. / ( 10. ** (2*(indice_epicont(j)-1)*ishelfal) )
        if (indice_epicont(j).ne.1.and.ishelfal.eq.1) then
            shelfal=0.01
        else
            shelfal=1.
        endif
        finorgC(j)=fbioC(j)*rC_Corg(j)*shells*shelfal
        finorgP(j)=finorgC(j)/1000.*0. ! pas de P ds les coquilles !

    enddo
    return
    end
