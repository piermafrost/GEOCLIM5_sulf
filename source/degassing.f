    subroutine degassing(t)
!   ***********************
    use constante, only: Gt_to_mol, tstart_deg, xhydLiin
    implicit none
    include 'combine_foam.inc'
    
    fvol=volin*clo    !*((-0.5/1d7)*t+1)


    do j0=1,ndeep
        j = jbox_deep(j0)
        fmor(j)=xMORin*clo/float(ndeep)
        FhydLi(j)=xhydLiin*clo/float(ndeep)
    end do
    do j0=1,nnodeep
        j = jbox_nodeep(j0)
        fmor(j)=0.
        FhydLi(j)=0.
    end do


    ftrap=0.
    !TRAP DEGASSING
    if (ipeak.eq.1) then
        do j=1,n_peaks    
            if (t.gt.tstart_deg+pulse_start(j).and.t.lt.tstart_deg+pulse_end(j)) then
                peak_duration(j)=pulse_end(j)-pulse_start(j)      !duration of one degassing peak
                ftrap=amount_peak(j)*Gt_to_mol/peak_duration(j)   !degassing in moles/yr
                dctrap=dc13_peak(j)
            endif
        end do
    endif
    

    fanthros=0.


    return
    end
