subroutine read_oceanic_temperature(IOunit, co2_axis, Toceclimber)

    implicit none
    include 'shape.inc'
    integer, intent(in) :: IOunit
    double precision, dimension(nclimber), intent(in) :: co2_axis
    double precision, dimension(nbasin, nclimber), intent(out) :: Toceclimber
    double precision, dimension(nclimber) :: SST, deepT
    double precision :: co2
    integer :: k
    character(len=10) :: line

    ! Try reading temperature generating mode (string) in input file
    read(unit=IOunit, fmt=*) line

    if (line == 'parametric') then
    
        ! Parametric CO2-temperature relationship
        ! ---------------------------------------

        !SST = 15. !6.8161*dlog(co2_axis*280.) - 22.089

        ! Pre-industrial
        SST = 20 + 4.88504*dlog(co2_axis)
        deepT = -6.58807 + 6.89026*dlog(co2_axis)
        ! Fixed oceanic temperature
        !SST = 20
        !deepT = -6.58807

        ! 265Ma
        !SST = 15.6279 + 4.88504*dlog(co2_axis)
        !deepT = -6.58807 + 6.89026*dlog(co2_axis)

        ! 365Ma
        !SST = 12.5071 + 5.98331*dlog(co2_axis)
        !deepT = -30.8490 + 16.9424*dlog(co2_axis)

        ! -----------------------------------------------------
        ! SURFACE:
        Toceclimber(1,:) = deepT  !SST-18
        Toceclimber(3,:) = SST
        Toceclimber(6,:) = SST
        Toceclimber(8,:) = deepT  !SST-18

        ! -----------------------------------------------------
        ! DEEP
        Toceclimber(2,:) = deepT
        Toceclimber(4,:) = deepT
        Toceclimber(5,:) = deepT
        Toceclimber(7,:) = Toceclimber(6,:)
        Toceclimber(9,:) = deepT


    else

      ! if input file does not contain the word "parametric",
      ! expect direct temperature reading, one CO2 level per line, in DECREASING order.

      rewind(IOunit)

      do k = nclimber,1,-1
          read(unit=IOunit, fmt=*) co2, Toceclimber(1:nbasin-1,k)
          !                   skip atmosphere box --^^^^^^^^^^
          if (abs(co2/280. - co2_axis(k)) > 1d-6) then
              print *
              print *, 'Error while reading oceanic temperature'
              print *, 'CO2 mismatch for level #', k
              write(*, fmt='(A,F10.2,A,F10.2)') 'Got: ', real(co2), ', expect: ', 280*co2_axis(k)
              stop
          end if
      end do


    end if


    ! Security: avoid too low oceanic temperature because it makes the chemistry routines crash
    where (Toceclimber(1:nbasin-1,:) <= 1)  Toceclimber(1:nbasin-1,:) = 1d0


    ! Temperature conversion: deg C => Kelvin
    Toceclimber(1:nbasin-1,:) = Toceclimber(1:nbasin-1,:) + 273.15


    close(unit=IOunit)

end subroutine
