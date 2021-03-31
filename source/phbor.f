    subroutine phbor(ct,alk,ah,sal,dens,ak1,ak2,akb)
!   ------------------------------------------------
    implicit none
    double precision ct,alk,ah,sal,dens,ak1,ak2,akb
    double precision cl,tbor,ba,sa,a,b,c,fah,dfah,residu,err

    cl=(sal-0.030)/1.8050
    tbor=2.19e-5*cl*dens
    ba=tbor/alk
    sa=ct/alk
    a=(akb*(1.-ba)+ak1*(1.-sa))*1.e+5
    b=(akb*(1.-ba-sa)+ak2*(1.-sa-sa))*ak1*1.e+10
    c=(1.-ba-sa-sa)*ak1*ak2*akb*1.e+15

!   methode de newton-raphson amelioree

    ah=max(abs(c),1.+abs(b),1.+abs(a))
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)

!   commencer par des pas doubles jusqua ce que fah<0

 1    continue
    ah=ah-2.*fah/dfah
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)
    if(fah.gt.0.)go to 1

!   continuer par newton-raphson classique

 2    continue
    residu=-fah/dfah
    ah=ah+residu
    fah=c+ah*(b+ah*(a+ah))
    dfah=b+ah*(a+a+3.*ah)
    err=abs(residu/ah)
    if(err.gt.1.e-3)go to 2

    ah=ah*1.e-5
    return
    end
