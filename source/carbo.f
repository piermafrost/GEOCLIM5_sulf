    subroutine carbo(Patm,temperature,runoff,xmass_carbonate,j)
!   ======================================================
    implicit none

    include 'combine_foam.inc'

    Pmax=1.+0.302483*runoff**0.8
    !Patm=var_diss(7,nbasin)/0.0508d+18 ! obsolete, Patm is now an input argument
    pco2soil=Patm+Pmax/(1+dexp(1.315-0.119*(temperature)))
    SO4soil=SO4_conc(j)

    temper=temperature+273.15

    x1=1.d-7
    x2=1.d-1
    xacc=1.d-10
!    if (value_keep1(j).ne.0) then
!        xn=value_keep1(j)
!    else
    xn=1.d-4
!    endif

!    call rtbis(x1,x2,xacc,out,temper)
    call newton(xn,fx,dfx,xacc,temper)
    
    out=xn
!    value_keep1(j)=xn

    h2co3soil=betasoil*pco2soil
    hco3soil=ak1soil*h2co3soil/out
    co3soil=hco3soil*ak2soil/out
    Calcium=0.5*(hco3soil+2*co3soil+2.*SO4soil)  !mol/m3
    xmass_carbonate=Calcium !mol/m3

    return
    end
