    subroutine chimie(co2t,alk,ch,h2co3,hco3,co3,ak1,ak2)
!   -----------------------------------------------------
    implicit none
      
    double precision co2t,alk,ch,h2co3,hco3,co3,ak1,ak2

    hco3=co2t/(ch/ak1+1.+ak2/ch)
    co3=ak2*hco3/ch
    h2co3=hco3*ch/ak1
    return
    end
