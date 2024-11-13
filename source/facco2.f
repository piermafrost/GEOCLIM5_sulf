    subroutine facco2(time,pco2,fco2,fe)
!   ------------------------------------
    implicit none
    include 'combine_foam.inc'
    
    pmin=0.2
    prmax=4.
    aco2=0.3
    pco2pal=var_diss(7,nbasin)
     
!   volk, 94
!   beta=1. !full vegetation ! pco2 en pal
!   phalf=(prmax-1.)*(1.-pmin)
!   prprma=(pco2pal-pmin)/(phalf+pco2pal-pmin)
!   prter=beta*prprma*prmax
!   peff=pco2pal/10.+0.9*prter
!   fco2v=peff**aco2
!   berner, 94
    fco2v=(2*pco2pal/(1+pco2pal))**0.4
!   gymnosperms value:
    fev=0.75
    fev=1.0

!   volk, 94 
!   no vegetation -> pression ds le sol, peff
!   peff=pco2pal
!   fco2nv=peff**aco2
!   berner, 94
    fco2nv=pco2pal**0.5 ! nv no vegetation
    fenv=0.15
    fracveg=fvegin
    fco2=(1.-fracveg)*fco2nv+fracveg*fco2v
    fco2=1
    if (berner.eq.1.) then
        fe=(1.-fracveg)*fenv+fracveg*fev
    else
        fe=1.
    endif
    return
    end
