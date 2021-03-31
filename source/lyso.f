    subroutine lyso(omc0,tem,sal,dens,dplyc,dplysa)
!   ------------------------------------------------
    implicit none
    double precision omc0,tem,sal,dens,dplyc,dplysa
    double precision grav,tcel,rt,cc,dvc,dkc,aa,bb,reali,pratmc
    double precision pratma

    grav=9.806
    tcel=tem-273.15
    rt=82.056*tem
    cc=-dlog(omc0)
    dvc=-(48.76-0.5304*tcel)
    dkc=-1.e-3*(11.76-0.3692*tcel)
    aa=0.5*dkc/rt
    bb=-dvc/rt
    reali=bb*bb-4.*aa*cc
    pratmc=(-bb+dsqrt(reali))/(2.*aa)
    dplyc=1.e-3*pratmc*1.013e5/(dens*grav)
    cc=-dlog(omc0/1.5)
    dvc=-(48.76-0.5304*tcel-2.8)
    bb=-dvc/rt
    reali=bb*bb-4.*aa*cc
    pratma=(-bb+dsqrt(reali))/(2.*aa)
    dplysa=1.e-3*pratma*1.013e5/(dens*grav)

    return
    end
