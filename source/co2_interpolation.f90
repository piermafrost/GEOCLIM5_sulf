module co2_interpolation
implicit none

contains

    subroutine climate_interpolation(co2lev)

        double precision, intent(in):: co2lev
        include 'combine_foam.inc'
        double precision:: xi
        integer:: nclimber_minus_1

        !------------------------------------------------------------------------------!
        ! If only one CO2 level => assume constant initial CO2 for silicate weathering !
        !------------------------------------------------------------------------------!
        if (nclimber==1) then

            ! Fixed CO2 mode => keep climate at level #1
            do j0=1,ncontpxl
                j = list_cont_pixel(j0)
                Tclim(j)   = Tairclimber(j,1)
                runclim(j) = Runclimber(j,1)
            end do
            !! "Fake" climate generator:
            !!   uniform 4°C warming and 25% runoff increase per doubling of CO2
            !do j0=1,ncontpxl
            !    j = list_cont_pixel(j0)
            !    Tclim(j)   = Tairclimber(j,1) + (4./log(2.)) * log(co2lev/co2climber(1))
            !    runclim(j) = Runclimber(j,1) * (1 + (0.25/log(2.)) * log(co2lev/co2climber(1)))
            !end do


        !--------------------------------------------------------------------!
        ! Otherwise, interpolate temperature and runoff at current CO2 level !
        !--------------------------------------------------------------------!
        else

            call find_co2_interval(co2climber, co2lev, k1, k2, xi)

            ! Case "out of CO2 range"
            if (co2lev < co2climber(1) .or. co2lev > co2climber(nclimber)) then

                do j0=1,ncontpxl
                    j = list_cont_pixel(j0)
                    ! "linear" extrapolation
                    Tclim(j) = (1-xi)*Tairclimber(j,k1) + xi*Tairclimber(j,k2)
                    runclim(j) = (1-xi)*Runclimber(j,k1) + xi*Runclimber(j,k2)
                    runclim(j) = max(  0d0  ,  (1-xi)*Runclimber(j,k1) + xi*Runclimber(j,k2)  )
                    ! security to avoid negative runoff in the extrapolation
                end do

            else ! Normal case

                do j0=1,ncontpxl
                    j = list_cont_pixel(j0)
                    ! "linear" interpolation
                    Tclim(j) = (1-xi)*Tairclimber(j,k1) + xi*Tairclimber(j,k2)
                    runclim(j) = (1-xi)*Runclimber(j,k1) + xi*Runclimber(j,k2)
                end do

            endif


            ! Spatialized interpolation coefficient A(x,y) and B(x,y), logarithic interpolation
            ! 
            !do j0=1,ncontpxl
            !    j = list_cont_pixel(j0)
            !    Tclim(j)   = ATe(j)+BTe(j)*dlog(co2lev)
            !    runclim(j) = ARu(j)+BRu(j)*Tclim(j)
            !    if (runclim(j).lt.0) runclim(j)=0
            !end do

        end if

    end subroutine


    !==============================================================================================!
    

    subroutine find_co2_interval(co2_axis, co2, kinf, ksup, xi)

        include 'shape.inc'
        double precision, dimension(nclimber), intent(in):: co2_axis
        double precision, intent(in)::                      co2
        integer, intent(out)::                              kinf, ksup
        double precision, intent(out)::                     xi
        integer::                                           k

        kinf = 1
        ksup = nclimber
        ! find by dichotomy:
        do while ( ksup-kinf > 1 )
            k = kinf + (ksup-kinf)/2
            if (co2 >= co2_axis(k)) then
                kinf = k
            else
                ksup = k
            end if
        end do

        ! interpolation coefficient => depends on interpolation mode (linear or logarithmic):
        xi = interp_coeff(co2, co2_axis(kinf), co2_axis(ksup))

    end subroutine


    !==============================================================================================!


    function interp_coeff(x, x0, x1)

        double precision:: interp_coeff
        double precision, intent(in):: x, x0, x1
        include 'coupler.inc'

        select case (interpolation_mode)
          case ('linear')
            interp_coeff = (x - x0) / (x1 - x0)
          case ('log')
            interp_coeff = log(x/x0) / log(x1/x0)
          case default
            print *
            print *, 'ERROR in "source/coupler.inc": illegal interpolation mode "'//interpolation_mode//'"'
            print *, 'Expect "linear" or "log"'
            stop
        end select

    end function
    


end module
