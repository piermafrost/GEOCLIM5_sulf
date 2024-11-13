    subroutine anoxic
!   ********************
    use constante, only: cp_oxic, cp_anoxic
!   C/P molar ratio: cp_oxic = 200 , cp_anoxic = 4000 (traditionnaly)
!   => ratio of what is buried in POC
    implicit none
    include 'combine_foam.inc'

    do i=1,nbasin
        call DOfA(var(11,i),DOA(i))
        cp_burial(i)=(cp_oxic*cp_anoxic)/ &
                   ((1-DOA(i))*cp_anoxic+DOA(i)*cp_oxic)
        if (cp_burial(i).gt.cp_anoxic) then
            cp_burial(i)=cp_anoxic
        endif
    enddo
    return
    end
