    subroutine org_dep(t)
!   *********************
    use constante, only: PIC_mmass, POC_mmass, rho_sed, gotoshelf, hml, hsr, betahml, gammahsr, xmethan!, fSO4cste
    implicit none
    include 'combine_foam.inc'


    ! Compute sedimentation flux (m3/yr) on each box
    ! ----------------------------------------------

    fexport = 0
    ! Sedimentation on continental shelf
    do j0 = 1,nsediepicontsurf
        j = jbox_sediepicontsurf(j0)
        fin0 = ( sedim_fract(j)*tss + &
                 gotoshelf*(fsink_inorg(j)*PIC_mmass*var(10,j) + fsink(j)*POC_mmass*var(9,j))*vol(j) ) / rho_sed
        ! Michaelis-like saturation: export what exceeds sedimentation capacity
        fin_sed(j) = fin0 / (1 + fin0/sedim_capacity(j))
        fexport = fexport + fin0-fin_sed(j)
    end do

    fsedim = fexport
    fexport = 0
    ! Sedimentation in epicontinental deep
    do j0 = 1,nsediepicontnosurf
        j = jbox_sediepicontnosurf(j0)
        fin0 = sedim_fract(j)*fsedim + &
               (fsink_inorg(j)*PIC_mmass*var(10,j)*vol(j) + fsink(j)*POC_mmass*var(9,j)*vol(j)) / rho_sed
        ! Michaelis-like saturation: export what exceeds sedimentation capacity
        fin_sed(j) = fin0 / (1 + fin0/sedim_capacity(j))
        fexport = fexport + fin0-fin_sed(j)
    end do

    fsedim = fexport
    ! Sedimentation in deep open ocean
    do j0 = 1,nsedinoepicont
        j = jbox_sedinoepicont(j0)
        fin_sed(j) = sedim_fract(j)*fsedim + &
                     (fsink_inorg(j)*PIC_mmass*var(10,j)*vol(j) + fsink(j)*POC_mmass*var(9,j)*vol(j)) / rho_sed
        ! Note: deep basins have infinite sedimentation capacity
    end do


!   Francois and Walker formalism, with TSS-dependent sedim rate:
!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do j0=1,nsedi
        j = jbox_sedi(j0)

!       sedimentation rate (m/yr)
        ws(j) = fin_sed(j) / surf_sedi(j)
!       input at the sediment top (mole(C)/a/m2):
        fin(j) = (1-indice_surface(j)*(1-gotoshelf)) * fsink(j)*var(9,j)*vol(j)/surf_sedi(j)
!       Corg at the basis of the bioturbated layer (mixed layer):
        Corg_hml(j) = fin(j)/(ws(j) + betahml*var(11,j)*hml)
        !! no oxygen feedback:
        !Corg_hml(j) = fin(j)/(ws(j) + betahml*0.23406*hml)
!       Corg at the basis of the sulfate reduction zone
        Corg_hsr(j) = (ws(j)*Corg_hml(j))/(ws(j) + gammahsr*var(20,j)*hsr) ! var(20,:) : [SO4^2-]

!       Convert in specific fluxes (per squared meters):
        fodc_noSR_m2 = ws(j)*Corg_hml(j)
        fodcm2(j)    = ws(j)*Corg_hsr(j)
        ! Sulfate-Reduction flux:
        fSulfRed_m2  = (1./2.)*(fodc_noSR_m2 - fodcm2(j))
        !              ^^^^^^^ 1 S reduced for 2 C oxidized

!       Oxygen flux: not affected by sulfate-reduction, except for the Fe^2+ part
        fO2_odc_m2   = fodc_noSR_m2 - (1./8.)*fSulfRed_m2
!                                   ^^^^^^^^^^^^^^^^^^^^^   for 1 S reduced, 1/2 Fe^2+ released => 1/8 O2 eventually consumed

!       Substract C loss through methanogenesis 
        fO2_odc_m2   = fO2_odc_m2 - (1-xmethan)*fodcm2(j) ! leaking methane will eventually be oxidized by O2
        fodcm2(j)    = xmethan*fodcm2(j)

!       Total fluxes:
        fodc(j)     = fodcm2(j)*surf_sedi(j)*clo
        fodp(j)     = fodc(j)/cp_burial(j)*phosss
        fSulfRed(j) = fSulfRed_m2*surf_sedi(j)*clo
        fO2_odc(j)  = fO2_odc_m2*surf_sedi(j)*clo

        ! Burial efficiency of C
        Corg_BE(j) = fodcm2(j)/fin(j)

    end do


    ! Null fluxes in basin not sedimenting
    do j0=1,nnosedi
        j = jbox_nosedi(j0)
        ws(j)       = 0
        fin_sed(j)  = 0
        fodc(j)     = 0
        fodp(j)     = 0
        fO2_odc(j)  = 0
        fSulfRed(j) = 0
        Corg_BE(j)  = 0
    end do


    return
    end
