module GCM_io_module
! contains subroutines handle to GCM outputs (ie, annual climatologies)
!   1. read the paths of configuration variables in the config file (config/GCM_io_condition)
!   2. read the continental variables (area, landfrac, temperature and runoff) from GCM outputs

    implicit none
    integer, parameter, private:: fmaxlen=500, vmaxlen=200
    double precision, parameter, private:: axis_relat_accuracy=1d-5
    ! maximum allowed relative difference between axis (lon, lat) from different files

    contains


    subroutine load_climatology(fID)

        use netcdf

        include 'combine_foam.inc'

        integer, intent(in):: fID
        character(len=fmaxlen), dimension(2):: area_file
        character(len=fmaxlen), dimension(nclimber):: climo_file, gmst_climo_file
        character(len=vmaxlen), dimension(2):: area_x_varname, area_y_varname, area_varname
        character(len=vmaxlen):: x_varname, y_varname, temp_varname, runf_varname, var_units, vname
        character(len=vmaxlen):: gmst_x_varname, gmst_y_varname, glob_temp_varname
        double precision:: fillvalue

        double precision, dimension(2+2*nclimber, nlon):: all_x
        double precision, dimension(2+2*nclimber, nlat):: all_y
        double precision, dimension(nlon, nlat):: totarea2D, landarea2D
        double precision, dimension(nclimber, nlon, nlat):: temperature2D, runoff2D, glob_temperature
        integer:: ierr, nax


        print *
        print *
        print *, 'Read GCM conditions file and load land inputs from GCM climatology files'
        print *, '------------------------------------------------------------------------'

        ! get netCDF file info
        call read_GCM_condition(fID, &
        !                       <><><><><>
                                co2climber, &
        !                       <><><><><>
                                area_file, area_x_varname, area_y_varname, area_varname, &
                                climo_file, x_varname, y_varname, temp_varname, runf_varname, &
                                gmst_climo_file, gmst_x_varname, gmst_y_varname, glob_temp_varname)


        ! get total area
        call load_variable('area', area_varname(1), area_x_varname(1), area_y_varname(1), &
                           single_input_file=area_file(1), varout2D=totarea2D, x=all_x(1,:), y=all_y(1,:), &
                           fill_missing=.true.)

        ! get land fraction/area
        call load_variable('landarea', area_varname(2), area_x_varname(2), area_y_varname(2), &
                           single_input_file=area_file(2), varout2D=landarea2D, x=all_x(2,:), y=all_y(2,:), &
                           totarea=totarea2D, fill_missing=.true.)

        ! get temperature
        call load_variable('temperature', temp_varname, x_varname, y_varname, &
                           multiple_input_file=climo_file, varout3D=temperature2D, &
                           xvec=all_x(3:2+nclimber,:), yvec=all_y(3:2+nclimber,:), &
                           landarea=landarea2D, fill_missing=.false., fillval_handling_option=ERROR_HANDLING_OPTION(2))
        ! note: ERROR_HANDLING_OPTION(2) tells the code what to do if it found missing-value where area > 0

        ! get runoff
        call load_variable('runoff', runf_varname, x_varname, y_varname, &
                           multiple_input_file=climo_file, varout3D=runoff2D, &
                           landarea=landarea2D, fill_missing=.false., fillval_handling_option=ERROR_HANDLING_OPTION(2))

        ! get optional global temperature
        if (glob_temp_varname /= '-') then
            call load_variable('temperature', glob_temp_varname, gmst_x_varname, gmst_y_varname, &
                            multiple_input_file=gmst_climo_file, varout3D=glob_temperature, &
                            xvec=all_x(3+nclimber:2+2*nclimber, :), yvec=all_y(3+nclimber:2+2*nclimber, :), &
                            fill_missing=.true.)
            ! number of axis to check
            nax = 2+2*nclimber
        else
            ! number of axis to check
            nax = 2+nclimber
        end if


        call check_axis(all_x(1:nax,:), all_y(1:nax,:), ERROR_HANDLING_OPTION(1))
        ! note: ERROR_HANDLING_OPTION(1) tells the code what to do if it find axis mismatch



        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!
        !       Put the loaded variables in GEOCLIM global variables       !
        !
        ref_x_axis = all_x(1,:)
        ref_y_axis = all_y(1,:)
        ! unravel 2D-arrays in 1D-array, the first dimension (x) begin the most rapidly varying => "natural" order
        areaEarth   = reshape(totarea2D,  shape=(/npixel/))
        areaclimber = reshape(landarea2D, shape=(/npixel/))
        !    Transpose CO2 axis, from 1st dimension to last dimension --vvvvvvvvvvvvv
        Tairclimber = reshape(temperature2D, shape=(/npixel,nclimber/), order=(/2,1/))
        Runclimber  = reshape(runoff2D,      shape=(/npixel,nclimber/), order=(/2,1/))
        ! Optional: Compute GMST
        if (glob_temp_varname /= '-') then
            totarea = sum(totarea2D)
            do k=1,nclimber
                GMSTclimber(k) = sum(glob_temperature(k,:,:)*totarea2D) / totarea
            end do
        else
            GMSTclimber = -1d99 ! "internal" missing-value -> will be replaced by "output" missing-value
        end if
        !                                                                  !
        !<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>!


    end subroutine


    !======================================================================!


    subroutine read_GCM_condition(fID, &
                                  CO2_levels, &
                                  area_file, area_x_varname, area_y_varname, area_varname, &
                                  climo_file, x_varname, y_varname, temp_varname, runf_varname, &
                                  gmst_climo_file, gmst_x_varname, gmst_y_varname, glob_temp_varname)

        use utils, only: add_path, read_comment

        integer, intent(in):: fID
        double precision, dimension(:), intent(out):: CO2_levels
        character(len=fmaxlen), dimension(2), intent(out):: area_file
        character(len=vmaxlen), dimension(2), intent(out):: area_x_varname, area_y_varname, area_varname
        character(len=fmaxlen), dimension(:), intent(out):: climo_file, gmst_climo_file
        character(len=vmaxlen), intent(out):: x_varname, y_varname, temp_varname, runf_varname
        character(len=vmaxlen), intent(out):: gmst_x_varname, gmst_y_varname, glob_temp_varname
        character(len=50):: tag

        integer:: nlevels
        integer:: k

        nlevels = size(CO2_levels)

        ! Get CO2 values
        call read_comment(fID)
        read(unit=fID, fmt=*) CO2_levels

        ! Get area and landfrac info
        do k = 1,2
            call read_comment(fID)
            read(unit=fID, fmt=*) area_file(k), area_x_varname(k), area_y_varname(k), area_varname(k)
            call add_path(area_file(k))
        end do

        ! Get climatology variables names
        call read_comment(fID)
        read(unit=fID, fmt='(A)') x_varname 
        call read_comment(fID)
        read(unit=fID, fmt='(A)') y_varname 
        call read_comment(fID)
        read(unit=fID, fmt='(A)') temp_varname 
        call read_comment(fID)
        read(unit=fID, fmt='(A)') runf_varname 

        ! Get climatology file paths
        call read_comment(fID)
        read(unit=fID, fmt=*) tag
        if (tag /= '<<--START-->>') then
            print *
            print *, 'Error: Expected "START" tag not found while reading GCM output file names'
            stop
        end if
        !
        do k = 1,nlevels
            call read_comment(fID)
            read(unit=fID, fmt=*) climo_file(k)
            call add_path(climo_file(k))
        end do
        !
        call read_comment(fID)
        read(unit=fID, fmt=*) tag
        if (tag /= '<<--STOP-->>') then
            print *
            print *, 'Error: Expected "STOP" tag not found while reading GCM output file names'
            print *, '       => number of input GCM files greater than size of CO2 axis'
            stop
        end if

        ! Get global temperature climatology variables names
        call read_comment(fID)
        read(unit=fID, fmt='(A)') gmst_x_varname 
        call read_comment(fID)
        read(unit=fID, fmt='(A)') gmst_y_varname 
        call read_comment(fID)
        read(unit=fID, fmt='(A)') glob_temp_varname 

        if (glob_temp_varname /= '-') then

            ! Get global temperature climatology file paths
            call read_comment(fID)
            read(unit=fID, fmt=*) tag
            if (tag /= '<<--START-->>') then
                print *
                print *, 'Error: Expected "START" tag not found while reading GCM output file names for GMST'
                stop
            end if
            call read_comment(fID)
            read(unit=fID, fmt=*) tag
            if (tag == '<<--STOP-->>') then ! No file specified => keep same as "main" climatology
                gmst_x_varname = x_varname
                gmst_y_varname = y_varname
                gmst_climo_file = climo_file
            else
                backspace(unit=fID)
                !
                do k = 1,nlevels
                    call read_comment(fID)
                    read(unit=fID, fmt=*) gmst_climo_file(k)
                    call add_path(gmst_climo_file(k))
                end do
                !
                call read_comment(fID)
                read(unit=fID, fmt=*) tag
                if (tag /= '<<--STOP-->>') then
                    print *
                    print *, 'Error: Expected "STOP" tag not found while reading GCM output file names for GMST'
                    print *, '       => number of input GCM files greater than size of CO2 axis'
                    stop
                end if
            end if

        end if

    end subroutine


    !======================================================================!


    subroutine load_variable(internal_varname, varname, x_varname, y_varname,         &
                             single_input_file, multiple_input_file,                  &
                             varout2D, varout3D, x, y, xvec, yvec, landarea, totarea, &
                             fill_missing, fillval_handling_option                    )
        use netcdf
        include 'shape.inc'

        character(len=*), intent(in):: internal_varname
        character(len=vmaxlen), intent(in):: varname, x_varname, y_varname
        character(len=fmaxlen), intent(in), optional:: single_input_file
        character(len=fmaxlen), dimension(nclimber), intent(in), optional:: multiple_input_file
        double precision, dimension(nlon,nlat), intent(inout), optional:: landarea
        double precision, dimension(nlon,nlat), intent(in), optional:: totarea
        double precision, dimension(nlon,nlat), intent(out), optional:: varout2D
        double precision, dimension(nclimber, nlon,nlat), intent(out), optional:: varout3D
        double precision, intent(out), optional:: x(nlon), y(nlat), xvec(nclimber,nlon), yvec(nclimber,nlat)
        integer, intent(in), optional:: fillval_handling_option
        logical, intent(in), optional:: fill_missing

        double precision, dimension(nlon,nlat):: dummyvar2D
        double precision, dimension(nclimber,nlon,nlat):: dummyvar3D
        double precision:: fillvalue
        character(len=vmaxlen):: var_units, vname
        character(len=1):: oper
        logical:: loc_fill_missing
        integer:: ierr, check, loc_fillval_handling_option
        integer:: k, n, noperations


        ! How to handle missing-value
        if (present(fill_missing)) then
            loc_fill_missing = fill_missing
        else
            loc_fill_missing = .false. ! default behaviour
        end if
        !
        ! Note: checking missing-points is done only if fill_missing==.false. and "landarea" argument is given
        if (present(fillval_handling_option)) then
            loc_fillval_handling_option = fillval_handling_option
        else
            loc_fillval_handling_option = -1 ! default value => ask user
        end if


        ! Check consistency between 2D variable or 3D (ie, multiple CO2 levels) variable
        ! => the right "input_file" and "var" arguments must be given.
        check = -11
        if (present(varout2D)) then
            check = check + 1
        end if
        if (present(single_input_file)) then
            check = check + 10
        end if
        if (present(varout3D)) then
            check = check + 2
        end if
        if (present(multiple_input_file)) then
            check = check + 9
        end if
        if (check /= 0) then
            print *
            print *, 'INTERNAL ERROR: inconsistent set of optional variable passed to the subroutine'
            print *, '"load_variable", in module "GCM_io_module"'
            stop
        end if


        ! Instead of reading just 1 variable name, they may be several variables,
        ! with arithmetic operation between them (addition or substraction)
        noperations = get_arithmetic_operations(varname)
        ! => write the scratch files:
        !        unit=334: list of individual variables
        !        unit=335: list of mathematics operators (+ or -)


        if (present(varout2D)) then

            ! Initialization
            varout2D = 0d0

            ! Loop to get variables and perform arithmetic operations:
            do n = 1,noperations

                ! load current variable
                read(unit=334, fmt=*) vname
                if (loc_fill_missing) then
                    ! load variable
                    if (present(x) .and. present(y)) then
                        call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                                x=x, y=y, varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                    else
                        call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                                varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                    end if
                    ! set var=0 on "missing" cells
                    if (ierr==NF90_NOERR) where (dummyvar2D==fillvalue) dummyvar2D = 0d0
                    ! check variable units
                    if (present(totarea)) then
                        call check_units(internal_varname, vname, dummyvar2D, var_units, totarea=totarea)
                    else
                        call check_units(internal_varname, vname, dummyvar2D, var_units)
                    end if
                else
                    ! load variable
                    if (present(x) .and. present(y)) then
                        call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                                x=x, y=y, varunits=var_units, fillval=fillvalue)
                    else
                        call load_netcdf_dble2D(single_input_file, x_varname, y_varname, vname, dummyvar2D, &
                                                varunits=var_units, fillval=fillvalue)
                    end if
                    ! check missingpoint
                    if (present(landarea)) then
                        call check_missingpoints(vname, single_input_file, dummyvar2D, fillvalue, landarea, &
                                                 loc_fillval_handling_option)
                    end if
                    ! check variable units
                    if (present(totarea)) then
                        call check_units(internal_varname, vname, dummyvar2D, var_units, fillvalue=fillvalue, totarea=totarea)
                    else
                        call check_units(internal_varname, vname, dummyvar2D, var_units, fillvalue=fillvalue)
                    end if
                end if

                ! perform arithmetic operation
                read(unit=335, fmt=*) oper
                select case (oper)
                    case ("+")
                        varout2D = varout2D + dummyvar2D
                    case ("-")
                        varout2D = varout2D - dummyvar2D
                end select

            end do


        else ! varout3D case => load var for each CO2 level

            ! Initialization
            varout3D = 0d0

            ! Loop to get variables and perform arithmetic operations:
            do n = 1,noperations

                ! load current variable
                read(unit=334, fmt=*) vname
                do k = 1,nclimber
                    if (loc_fill_missing) then
                        ! load variable
                        if (present(xvec) .and. present(yvec)) then
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(k,:,:), &
                                               x=xvec(k,:), y=yvec(k,:), varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                        else
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(k,:,:), &
                                                    varunits=var_units, fillval=fillvalue, fillval_iostat=ierr)
                        end if
                        ! set var=0 on "missing" cells
                        if (ierr==NF90_NOERR) where (dummyvar3D(k,:,:)==fillvalue) dummyvar3D(k,:,:) = 0d0
                        ! check variable units
                        if (present(totarea)) then
                            call check_units(internal_varname, vname, dummyvar3D(k,:,:), var_units, totarea=totarea)
                        else
                            call check_units(internal_varname, vname, dummyvar3D(k,:,:), var_units)
                        end if
                    else
                        ! load variable
                        if (present(xvec) .and. present(yvec)) then
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(k,:,:), &
                                                    x=xvec(k,:), y=yvec(k,:), varunits=var_units, fillval=fillvalue)
                        else
                            call load_netcdf_dble2D(multiple_input_file(k), x_varname, y_varname, vname, dummyvar3D(k,:,:), &
                                                    varunits=var_units, fillval=fillvalue)
                        end if
                        ! check missingpoint
                        if (present(landarea)) then
                            call check_missingpoints(vname, multiple_input_file(k), dummyvar3D(k,:,:), fillvalue, landarea, &
                                                     loc_fillval_handling_option)
                        end if
                        ! check variable units
                        if (present(totarea)) then
                            call check_units(internal_varname, vname, dummyvar3D(k,:,:), var_units, fillvalue=fillvalue, &
                                             totarea=totarea)
                        else
                            call check_units(internal_varname, vname, dummyvar3D(k,:,:), var_units, fillvalue=fillvalue)
                        end if
                    end if
                end do

                ! perform arithmetic operation
                read(unit=335, fmt=*) oper
                select case (oper)
                    case ("+")
                        varout3D = varout3D + dummyvar3D
                    case ("-")
                        varout3D = varout3D - dummyvar3D
                end select

            end do


        end if


        ! close scratch files
        close(unit=334)
        close(unit=335)


    end subroutine


    !======================================================================!


    function get_arithmetic_operations(varstring)
    ! => create scratch files 334 and 335 containing (respectively) the lists of variables and operators
    ! Note: the first operators (corresponding to the first variable) is automatically '+'
    ! return the number of variables/operations
        character(len=*), intent(in):: varstring
        integer:: get_arithmetic_operations
        character(len=1), dimension(2), parameter:: operators = (/'+', '-'/) ! Legal arithmetic operators
        integer:: i0, i1
        open(unit=334, status='scratch') ! individual variables
        open(unit=335, status='scratch') ! operators
        i0 = 1
        i1 = 1
        get_arithmetic_operations = 1
        write(unit=335, fmt='(A)') '+'
        do i1 = 1,len(varstring)
            if (any(operators == varstring(i1:i1))) then
                write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:i1-1)))
                write(unit=335, fmt='(A)') varstring(i1:i1)
                i0 = i1+1
                get_arithmetic_operations = get_arithmetic_operations + 1
            end if
        end do
        write(unit=334, fmt='(A)') trim(adjustl(varstring(i0:)))
        rewind(unit=334)
        rewind(unit=335)
    end function


    !======================================================================!


    subroutine load_netcdf_dble2D(fname, x_varname, y_varname, varname, var, x, y, varunits, fillval, fillval_iostat)

        use netcdf

        character(len=*), intent(in):: fname, x_varname, y_varname, varname
        double precision, dimension(:,:), intent(out):: var
        double precision, dimension(:), intent(out), optional:: x, y
        character(len=*), intent(out), optional:: varunits
        double precision, intent(out), optional:: fillval
        integer, intent(out), optional:: fillval_iostat
        integer:: nx, ny, ndim, k, n
        integer:: ierr, fid, varid, xdimid, ydimid, shp_ix, shp_iy
        integer, dimension(:), allocatable:: dimids, shp
        logical:: transp
        double precision, dimension(:,:), allocatable:: loc_var

        nx = size(var, 1)
        ny = size(var, 2)
        ! must also match the size of "x" and "y"

        ! open netCDF file
        ierr = nf90_open(fname, NF90_NOWRITE, fid)
        call nf90_check(ierr, 'Error while openning file '//trim(fname))

        ! Get axis data
        if (present(x)) then
            call load_axis(fname, fid, x_varname, xdimid, x)
        else
            call load_axis(fname, fid, x_varname, xdimid)
        end if
        if (present(y)) then
            call load_axis(fname, fid, y_varname, ydimid, y)
        else
            call load_axis(fname, fid, y_varname, ydimid)
        end if

        ! Get variable ID
        ierr = nf90_inq_varid(fid, varname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(varname)//'" ID in file "'//trim(fname)//'"')

        ! Get variable number of dimension
        ierr = nf90_inquire_variable(fid, varid, ndims=ndim)
        call nf90_check(ierr, 'Error while inquiring number of dim. of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        allocate(dimids(ndim))
        allocate(shp(ndim))

        ! Get variable dimension IDs
        ierr = nf90_inquire_variable(fid, varid, dimids=dimids)
        call nf90_check(ierr, 'Error while inquiring dimensions IDs of variable "'//trim(varname)//'" in file "'//trim(fname)//'"')

        ! Get variable shape (inquire length of all dimensions)
        do k = 1,ndim
            ierr = nf90_inquire_dimension(fid, dimids(k), len=shp(k))
            call nf90_check(ierr, 'Error while inquiring length of dimension in file "'//trim(fname)//'"')
        end do

        ! Check variable shape (must have exactly 2 non degenerated dimensions)
        if (count(shp > 1) /= 2) then
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
                     '" must have exactly 2 non-degenerated (ie, size-1) dimensions'
            stop
        else
            k = 1
            do while (shp(k)==1)
                k = k + 1
            end do
            shp_ix = k
            k = k + 1
            do while (shp(k)==1)
                k = k + 1
            end do
            shp_iy = k
        end if

        ! Check that variable is defined on the given dimensions
        if (dimids(shp_ix)==xdimid .and. dimids(shp_iy)==ydimid) then
            transp = .false.
            allocate(loc_var(nx,ny))
        elseif (dimids(shp_ix)==ydimid .and. dimids(shp_iy)==xdimid) then
            transp = .true.
            allocate(loc_var(ny,nx))
        else
            print *
            print *
            print *, 'Error: variable "'//trim(varname)//'" of file "'//trim(fname)// &
            '" is not defined on the given dimensions "'//trim(x_varname)//'" and "'//trim(y_varname)//'"'
            stop
        end if

        ! load variable
        if (ndim>2) then
            print *
            print *,                 'Variable "'//trim(varname)//'" of file "'//trim(fname)//'"'
            write(*, fmt='(A,I2,A)') '    Ignore ', ndim-2, ' degenerated (size-1) dimension(s).'
        end if
        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
        ierr = nf90_get_var(fid, varid, loc_var)
        call nf90_check(ierr, 'Error while getting variable "'//trim(varname)//'" of file "'//trim(fname)//'"')
        ! + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + !
        if (transp) then
            var = transpose(loc_var)
        else
            var = loc_var
        end if

        ! Get variable units (if asked)
        if (present(varunits)) then
            varunits = ''
            ierr = nf90_get_att(fid, varid, 'units', varunits)
            call nf90_check(ierr, 'Warning: unable to get attribute "units" of variable "'//trim(varname)//'" in file "' &
                                   //trim(fname)//'". Variable assumed to be dimensionless.', kill=.false.)
        end if

        ! Get variale fill-value (if asked)
        if (present(fillval)) then
            ierr = nf90_get_att(fid, varid, '_FillValue', fillval)
            if (present(fillval_iostat)) then
                fillval_iostat = ierr
            else
                call nf90_check(ierr, 'Error while getting attribute "_FillValue" of variable "'//trim(varname)//'" in file "' &
                                                                                                          //trim(fname)//'"')
            end if
        end if

        ! close netCDF file
        ierr = nf90_close(fid)
        call nf90_check(ierr, 'Error while closing file '//trim(fname), kill=.false.)

        deallocate(dimids)
        deallocate(shp)
        deallocate(loc_var)

    end subroutine


    !======================================================================!


    subroutine load_axis(fname, fid, axname, axdimid, ax)
        use netcdf
        character(len=*), intent(in):: fname, axname
        integer, intent(in):: fid
        integer, intent(out):: axdimid
        double precision, dimension(:), intent(out), optional:: ax
        integer:: varid, ierr

        ! Inquire dim ID
        ierr = nf90_inq_dimid(fid, axname, axdimid)
        call nf90_check(ierr, 'Error while inquiring dimension "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Inquire var ID
        ierr = nf90_inq_varid(fid, axname, varid)
        call nf90_check(ierr, 'Error while inquiring variable "'//trim(axname)//'" ID in file "'//trim(fname)//'"')
        ! Get variable
        if (present(ax)) then
            ierr = nf90_get_var(fid, varid, ax)
            call nf90_check(ierr, 'Error while getting variable "'//trim(axname)//'" in file "'//trim(fname)//'"')
        end if

    end subroutine


    !======================================================================!


    subroutine compare_units(units_string, known_units, conversion_factor, conversion_offset, passed)

        use physical_units, only: units

        character(len=*), intent(in):: units_string
        type(units), intent(in):: known_units
        double precision, intent(out):: conversion_factor, conversion_offset
        logical, intent(out):: passed
        integer:: k

        conversion_factor = 1d0
        conversion_offset = 0d0

        passed = .false.

        if (units_string == known_units%reference) then
            passed = .true.
        else
            do k = 1,known_units%naccepted
                if (units_string == known_units%accepted(k)%string) then
                    passed = .true.
                    conversion_factor = known_units%accepted(k)%conversion(1)
                    conversion_offset = known_units%accepted(k)%conversion(2)
                end if
            end do
        end if

    end subroutine

    !======================================================================!

    subroutine check_units(which_variable, varname, var, varunits, fillvalue, totarea, error_handling)

        use physical_units, only: units, area_units, fraction_units, temperature_units, runoff_units

        character(len=*), intent(in):: which_variable, varname, varunits
        double precision, dimension(:,:), intent(inout):: var
        double precision, dimension(:,:), intent(in), optional:: totarea
        double precision, intent(in), optional:: fillvalue
        integer, intent(in), optional:: error_handling
        type(units):: known_units
        double precision:: factor, offset
        logical:: multiply_by_area, passed, loop
        integer:: k, loc_error_handling, answer


        if (present(error_handling)) then
            loc_error_handling = error_handling
        else
            loc_error_handling = -1 ! default behaviour: ask user
        end if


        multiply_by_area = .false.

        ! case-dependent statement
        select case (which_variable)

            case ("area")
                known_units = area_units()

            case ("landarea")
                known_units = area_units()

            case ("temperature")
                known_units = temperature_units()

            case ("runoff")
                known_units = runoff_units()

            case default
                print *
                print *, 'INTERNAL ERROR in function "check_units" of module "GCM_io_module": unkown variable case "' &
                                                                                            //trim(which_variable)//'"'
                stop

        end select


        ! Compare given units to reference units, and get conversion values if not equals
        call compare_units(varunits, known_units, factor, offset, passed)
        ! in the case "landarea", also try "fraction" units
        if (which_variable == "landarea" .and. (.not. passed)) then
            call compare_units(varunits, fraction_units(), factor, offset, passed)
            multiply_by_area = passed
        end if



        ! Print message
        if (passed) then
            if (factor/=1 .or. offset/=0 .or. multiply_by_area) then
                print *
                print *, 'Automatic conversion of variable "'//trim(varname)//'":'
                print *, '    "'//trim(varunits)//'" => "'//trim(known_units%reference)//'"'
            end if
        else
            print *
            print *, 'WARNING: unkown units of variable "'//trim(varname)//'".'
            print *, '    got units:      "'//trim(varunits)//'"'
            print *, '    expected units: "'//trim(known_units%reference)//'"'
            ! Error handling:
            select case (loc_error_handling)
                case (-1) ! ask user interactively
                    print *
                    loop = .true.
                    do while (loop)
                        loop = .false.
                        print *, 'Enter one of the following options:'
                        print *, '    0: abort the program'
                        print *, '    2: ignore the issue and continue the execution'
                        read *, answer
                        select case (answer)
                            case (0); stop
                            case (2) ! do nothing
                            case default; loop=.true.
                        end select
                    end do

                case (0) ! abort the program
                    stop

                case (2) ! do nothing

                case default
                    print *
                    print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                    stop
            end select
        end if


        ! Conversion
        if (present(fillvalue)) then
            where (var/=fillvalue)  var = factor*(var + offset)
        else
            var = factor*(var + offset)
        end if

        ! Multiply value by cell area
        if (multiply_by_area) then
            if (present(totarea)) then
                print *, 'Variable "'//trim(varname)//'" multiplied by total area'
                if (present(fillvalue)) then
                    where (var/=fillvalue) var = var * totarea
                else
                    var = var * totarea
                end if
            else
                print *
                print *, 'INTERNAL ERROR: total area not passed to the function "check_units" in module "GCM_io_module".'
                print *, 'You have the right to be mad at the developer.'
                stop
            end if
        end if


    end subroutine


    !======================================================================!


    subroutine check_missingpoints(varname, filename, var, fillvalue, landarea, error_handling)

        character(len=*), intent(in):: varname, filename
        double precision, dimension(:,:), intent(in):: var
        double precision, intent(in):: fillvalue
        double precision, dimension(:,:), intent(inout):: landarea
        integer, intent(in):: error_handling
        logical, dimension(:,:), allocatable:: errormask
        double precision:: area_err, tot_landarea
        integer:: nerr, answer
        logical:: loop

        tot_landarea = sum(landarea)

        allocate(errormask(size(landarea,1), size(landarea,2)))

        errormask = (landarea>0 .and. var==fillvalue)

        nerr = count(errormask)
        area_err = sum(landarea, mask=errormask)

        ! If error found:
        if (nerr > 0) then

            print *
            write(*,'(A)')         ' WARNING: found missing-value on continental points (land area > 0) in variable "' &
                                    //trim(varname)//'" of file "'//trim(filename)//'"'
            write(*,'(A,I0,A,I0)') '     Number of problematic continent cells:  ', nerr, ' / ', count(landarea>0)
            write(*,'(A,E14.7)')   '     Total area of those cells (km2):        ', 1d-6*area_err
            write(*,'(A,E14.7)')   '     Which is a fraction of total land area: ', area_err/tot_landarea

            ! Error handling:
            select case (error_handling)
                case (-1) ! ask user interactively
                    print *
                    loop = .true.
                    do while (loop)
                        loop = .false.
                        print *, 'Enter one of the following options:'
                        print *, '    0: abort the program'
                        print *, '    1: remove erratic points (set area=0)'
                        print *, '    2: ignore the issue and continue the execution'
                        read *, answer
                        select case (answer)
                            case (0); stop
                            case (1); where (errormask) landarea=0
                            case (2) ! do nothing
                            case default; loop=.true.
                        end select
                    end do

                case (0) ! abort the program
                    stop

                case (1) ! automatic correction
                    print *
                    print *, 'Automatic correction => remove erratic points (set area=0)'
                    where (errormask) landarea=0

                case (2) ! do nothing

                case default
                    print *
                    print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                    stop

            end select

        end if

        deallocate(errormask)

    end subroutine


    !======================================================================!


    subroutine raise_axis_error(which_axis, nerr, axis_len, max_mismatch, error_handling)
        character(len=1), intent(in):: which_axis
        integer, intent(in):: nerr, axis_len
        double precision, intent(in):: max_mismatch
        integer, intent(in):: error_handling
        integer:: answer
        logical:: loop

        print *
        write(*,'(A)')         ' WARNING: found mismatch of '//which_axis//' axis between GCM input files.'
        write(*,'(A,I0,A,I0)') '     Number of points concerned / axis length:  ', nerr, ' / ', axis_len
        write(*,'(A,E14.7)')   '     Maximum mismatch found (in axis units):    ', max_mismatch

        ! Error handling:
        select case (error_handling)
            case (-1) ! ask user interactively
                print *
                loop = .true.
                do while (loop)
                    loop = .false.
                    print *, 'Enter one of the following options:'
                    print *, '    0: abort the program'
                    print *, '    2: ignore the issue and continue the execution'
                    read *, answer
                    select case (answer)
                        case (0); stop
                        case (2) ! do nothing
                        case default; loop=.true.
                    end select
                end do

            case (0) ! abort the program
                stop

            case (1) ! asked for automatic correction
                print *
                print *, 'Bad error handling option "1": cannot remove points with axis mismatch.'
                stop

            case (2) ! do nothing

            case default
                print *
                print *, 'INTERNAL ERROR: illegal error handling option:', error_handling
                stop

        end select

    end subroutine

    !======================================================================!

    subroutine check_axis(all_x, all_y, error_handling)
        double precision, dimension(:,:), intent(in):: all_x, all_y ! dim #1: files (temp, runoff, ...), dim #2: axis
        integer, intent(in):: error_handling
        double precision:: daxis, max_mismatch
        logical, dimension(:), allocatable:: errormask
        integer:: n, k, axlen, nerr

        n = size(all_x, 1) ! also = size(all_y, 1)


        ! x axis
        ! ------

        axlen = size(all_x, 2)
        allocate(errormask(axlen))
        errormask = .false.

        daxis = abs(all_x(1,2) - all_x(1,1))
        max_mismatch = 0d0

        do k = 2,n
            errormask = (errormask  .or.  maxval(abs(all_x(k,:) - all_x(1,:)))/daxis > axis_relat_accuracy)
            max_mismatch = max(max_mismatch, maxval(abs(all_x(k,:) - all_x(1,:))))
        end do
        nerr = count(errormask)
        if (nerr > 0) then
            call raise_axis_error('x', nerr, axlen, max_mismatch, error_handling)
        end if

        deallocate(errormask)


        ! y axis
        ! ------

        axlen = size(all_y, 2)
        allocate(errormask(axlen))
        errormask = .false.

        daxis = abs(all_y(1,2) - all_y(1,1))
        max_mismatch = 0d0

        do k = 2,n
            errormask = (errormask  .or.  maxval(abs(all_y(k,:) - all_y(1,:)))/daxis > axis_relat_accuracy)
            max_mismatch = max(max_mismatch, maxval(abs(all_y(k,:) - all_y(1,:))))
        end do
        nerr = count(errormask)
        if (nerr > 0) then
            call raise_axis_error('y', nerr, axlen, max_mismatch, error_handling)
        end if

        deallocate(errormask)

    end subroutine


    !======================================================================!


    subroutine nf90_check(ierr, message, kill)

        use netcdf
        integer, intent(in):: ierr
        character(len=*), intent(in):: message
        logical, optional, intent(in):: kill
        logical:: loc_kill

        if (present(kill)) then
            loc_kill = kill
        else
            loc_kill = .true.
        end if

        if (ierr/=NF90_NOERR) then
            print *
            print *, message
            print *, nf90_strerror(ierr)
            if (loc_kill) stop
        end if

    end subroutine


end module
