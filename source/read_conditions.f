    subroutine read_conditions
!   **************************
    implicit none
    include 'combine_foam.inc'


    read(3,*) volin
    read(3,*) xMORin
    read(3,*) fSO4_atmos
    read(3,*) fSO4_ocean
    read(3,*) fSO4_deep
    read(3,*) fCO2_atmos
    read(3,*) fCO2_deep
    read(3,*) ipeak
    read(3,*) i_rnd
    read(3,*) fauna
    read(3,*) iextinct
    read(3,*) extrinsic_flag
    read(3,*) iweb
    read(3,*) clo
    read(3,*) shells
    read(3,*) ishelfal
    read(3,*) phosss
    read(3,*) xnoorg
    read(3,*) temcondi
    read(3,*) ts
    read(3,*) tsta
    read(3,*) tfi
    read(3,*) xjump
    read(3,*) ijump_geogprint
    read(3,*) ageYprint
    read(3,*) oxy_acc_fact
    read(3,*) sulf_acc_fact
    read(3,*) isolver
    read(3,*) ijump_cont_weath
    read(3,*) ijump_DS_integration
    read(3,*) ijump_DS_print
    read(3,*) coupling_veget
    read(3,*) ijump_veget
    read(3,*) convert2ascii

    ! generate restart at end of run if asked (0 or negative time)
    if (ageYprint <= 0) ageYprint = tfi + 100d6

    return
    end
