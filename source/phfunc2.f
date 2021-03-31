    subroutine phfunc2(xv,fph,dfph)
!   --------------------------------
!   This routine calculate pH for a solution at equilibrium with
!   soil PCO2. Used to estimate phosphorus weathering

    implicit none
    include 'combine_foam.inc'
    
    pr=0.
    
    fph=-xv*xv*xv+1.d-8*xv+2*ak1soil*ak2soil*betasoil*pco2soil+ &
        ak1soil*xv*betasoil*pco2soil
    dfph=-3.*xv*xv+1.d-8+ak1soil*xv*betasoil*pco2soil
    
    return
    end
