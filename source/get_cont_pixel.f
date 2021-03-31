subroutine get_cont_pixel()
    implicit none
    include 'combine_foam.inc'

    j = 0
    do j0 = 1,npixel
        if (areaclimber(j0)>0) then
            j = j+1
            list_cont_pixel(j) = j0
        end if
    end do

    ncontpxl = j

end subroutine
