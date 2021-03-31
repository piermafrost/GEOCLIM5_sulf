    subroutine eqcte(tem,sal,pr,bc,ak1,ak2,akb,akc)
!   -----------------------------------------------
    implicit none
    double precision tem,sal,pr,bc,ak1,ak2,akb,akc,ak0
    double precision pratm,tcel,rt,dv1,dk1,arg,dv2,dk2,dvb,dkb
!   this subroutine returns the equilibrium constants in mol/m3
!   for the dissociation of h2co3 (ak1), hco3 (ak2) and h3bo3 (akb)
!   as a function of temperature (tem) in kelvin, salinity (sal)
!   in per mil and pressure (pr) in bar. it also returns
!   henrys law constant bc for co2,(weiss (1974) for co2 solubility
!   in water, in zeebe & wolf-gladrow -co2 in seawater: equilibrium,
!   kinetics, isotopes-) in mol/m3/pal, and the
!   saturation product (akc) for calcite at pr=0 in mol/m3*mol/m3.

    ak0=-60.2409+9345.17/tem+23.3585*dlog(tem/100) &
        +sal*(0.023517-0.00023656*tem+0.0047036*(tem*tem/10000))
    ak0=dexp(ak0)*1.e+3
    bc=280.e-6*ak0
    ak1=290.9097-14554.21/tem-45.0575*dlog(tem) &
        +dsqrt(sal)*(0.0221+34.02/tem)
    ak1=dexp(ak1)*1.e+3
    ak2=207.6548-11843.79/tem-33.6485*dlog(tem) &
        +dsqrt(sal)*(0.9805-92.65/tem)-0.03294*sal
    ak2=dexp(ak2)*1.e+3
    akb=148.0248-8966.90/tem-24.4344*dlog(tem) &
        +dsqrt(sal)*(0.5998-75.25/tem)-0.01767*sal
    akb=dexp(akb)*1.e+3
    akc=303.1308-13348.09/tem-48.7537*dlog(tem) &
        +dsqrt(sal)*(1.6233-118.64/tem)-0.06999*sal
    akc=dexp(akc)*1.e+6
    if(pr.ge.1.)then
        pratm=pr/1.013
        tcel=tem-273.15
        rt=82.056*tem
        dv1=-(25.50+0.151*(sal-34.8)-0.1271*tcel)
        dk1=-1.e-3*(3.08+0.578*(sal-34.8)-0.0877*tcel)
        arg=-dv1*pratm/rt+0.5*dk1*pratm*pratm/rt
        ak1=ak1*dexp(arg)
        dv2=-(15.82-0.321*(sal-34.8)+0.0219*tcel)
        dk2=-1.e-3*(-1.13+0.314*(sal-34.8)+0.1475*tcel)
        arg=-dv2*pratm/rt+0.5*dk2*pratm*pratm/rt
        ak2=ak2*dexp(arg)
        dvb=-(29.48-0.295*(sal-34.8)-0.1622*tcel+0.002608*tcel*tcel)
        dkb=-1.e-3*(2.84-0.354*(sal-34.8))
        arg=-dvb*pratm/rt+0.5*dkb*pratm*pratm/rt
        akb=akb*dexp(arg)
    endif
    return
    end
