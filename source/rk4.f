    subroutine rk4(y,dydx,n,x,h,yout)
!   *********************************
    implicit none
    integer, parameter :: nmax=4000
    integer n,i
    double precision y(n),dydx(n),yout(n),yt(n),dyt(n),dym(n)
    double precision h,hh,h6,x,xh

    hh=h*0.5
    h6=h/6.
    xh=x+hh

    call derivs(x,y,n,dydx)

    do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
   11 continue

    call derivs(xh,yt,n,dyt)

    do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
   12 continue

    call derivs(xh,yt,n,dym)

    do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
   13 continue

    call derivs(x+h,yt,n,dyt)

    do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
   14 continue

    return
    end
