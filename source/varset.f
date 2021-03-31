    subroutine varset(y,nv)
!   ***********************
    implicit none
    integer nv

    include 'combine_foam.inc'
!   **********************
!   DIC(1),Alk(2),PO4(3),Ca(4),Sr(5),SrPIC(6)
!   POP(7),PIP(8),POC(9),PIC(10),O2(11),PCO2(12)
!   DIC13C(13),PIC13C(14),POC13C(15),PCO2_13C(16)
!   SrIso(17),SrPICISO(18),7Li/6Li(19),Li(20)
!   **********************
    
    k=0
    
    !========
    do i=1,10
    !--------
        do j=1,nbasin
            k=k+1
            var(i,j)=y(k)/vol(k)
        end do
    end do
    
    !=====
    i = 11
    !-----
    do j0=1,nnosurface ! do the calculation only in no-surface boxes
        j = jbox_nosurface(j0)
        var(i,j) = y(k+j)/vol(k+j)
    end do
    k = k+nbasin
    
    !=====
    i = 12
    !-----
    do j=1,nbasin
        k=k+1
        var(i,j)=y(k)/vol(k)
    end do
    
    !=========
    do i=13,19
    !---------
        do j=1,nbasin
            k=k+1
            var(i,j)=y(k)
        end do
    end do
    
    !=====
    i = 20
    !-----
    do j=1,nbasin
        k=k+1
        var(i,j)=y(k)/vol(k)
    end do

    return
    end
