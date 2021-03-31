    subroutine derivs(x,y,nv,dydx)
!   ******************************
    implicit none

    integer nv
    include 'combine_foam.inc'


    call varset(y,nv)

!   Initialise les derives

    do i=1,nvar_real
        do j=1,nbasin
            diff(i,j)=0.
        enddo
    enddo
    

!    hard code
!    if (x.ge.1.d+4.and.x.le.3d+4) then
!        sluggish=0.5
!    else
        sluggish=1.
!    endif

!   Calcul

!   cccccccccccccccccccccccccccccccccccccccc
!   dissovled variables: i= 1-5 , 12 & 20
!   cccccccccccccccccccccccccccccccccccccccc
!   --------
    do i=1,5
        do j=1,nbasin
            do k=1,nbasin
                diff(i,j)=diff(i,j)+F(k,j)*var(i,k)*sluggish  &
                                   -F(j,k)*var(i,j)*sluggish
            enddo
            diff(i,j)=diff(i,j) + R(i,j)
        enddo
    end do
!   ----
    i=12
    do j=1,nbasin
        do k=1,nbasin
            diff(i,j)=diff(i,j)+F(k,j)*var(i,k)*sluggish  &
                               -F(j,k)*var(i,j)*sluggish
        enddo
        diff(i,j)=diff(i,j) + R(i,j)
    enddo
!   ----
    i=20  ! SO4^2-
    do j=1,nbasin
        do k=1,nbasin
            diff(i,j)=diff(i,j)+F(k,j)*var(i,k)*sluggish  &
                               -F(j,k)*var(i,j)*sluggish
        enddo
        diff(i,j)=diff(i,j) + R(i,j)
    end do

!       cccccccccccccccccccccccccccccccccccccccc
    i=11 ! Oxygen
!       cccccccccccccccccccccccccccccccccccccccc
    do j0=1,nnosurface
        j = jbox_nosurface(j0)
        do k=1,nbasin
            diff(i,j)=diff(i,j)+F(k,j)*var(i,k)*sluggish  &
                               -F(j,k)*var(i,j)*sluggish
        enddo
        diff(i,j)=diff(i,j) + R(i,j)
    end do
    do j0=1,nsurface !!! oxygN de surf. equilibr (diff(i,j)=0)
        j = jbox_surface(j0)
        diff(i,j)=R(i,j)
    end do
!       cccccccccccccccccccccccccccccccccccccccc
    do i=6,10 ! particulate variables
!       cccccccccccccccccccccccccccccccccccccccc
        do j=1,nbasin
            diff(i,j)=diff(i,j) + R(i,j)
        enddo
    end do

!       cccccccccccccccccccccccccccccccccccccccc
    do i=13,18 ! isotopic ratios
!       cccccccccccccccccccccccccccccccccccccccc
        do j=1,nbasin ! dissolved isotope are mixed in creades
            diff(i,j)=diff(i,j) + R(i,j)
        end do
    end do

!   Deplace les derivs ds un vector

    diff(11,10)=diff(11,10)*oxy_acc_fact !accelerating atm O2 (if oxy_acc_fact>1)
    diff(20,:) = sulf_acc_fact*diff(20,:) ! accelerating sulfate cycle (in all oceanic boxes)

    k=0
    do i=1,nvar_real
        do j=1,nbasin
            k=k+1
            dydx(k)=diff(i,j)
        enddo
    enddo

    return
    end
