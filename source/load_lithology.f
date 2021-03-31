subroutine load_lithology(ID)
    use netcdf
    use io_module, only: nf90check, check_landcells, check_axis, check_invalid
    use utils, only: add_path
    implicit none

    include 'combine_foam.inc'
    integer, intent(in):: ID
    integer:: fileid, varid, ierr
    integer, dimension(3):: dimids, dimlen
    character(len=50), dimension(3):: dimname
    double precision:: fillval
    logical:: got_fillval
    double precision, dimension(nlon):: loc_x_axis
    double precision, dimension(nlat):: loc_y_axis
    double precision, dimension(nlitho):: singlepixel_lithofrac


    print *
    print *
    print *, 'Load lithology map'
    print *, '------------------'
    print *

    !-----------!
    ! Load info !
    !-----------!
    rewind(unit=ID)

    ! Try direct litho frac reading (=> uniform lithology):
    read(unit=ID, fmt=*, iostat=ierr) singlepixel_lithofrac

    if (ierr==0) then

        ! print report:
        print *, '    Uniform lithology'
        print *

        do k = 1,nlitho
            litho_frac(k,:) = singlepixel_lithofrac(k)
        end do


    else ! expect netCDF file info

        rewind(unit=ID)
        read(unit=ID, fmt=*) dummychar, filename, varname, fillvalname
        call add_path(filename)


        !----------------!
        ! File openning: !
        !----------------!

        ierr = nf90_open(filename, NF90_NOWRITE, fileid)
        call nf90check(ierr,'Error while openning input file '//trim(filename))


        !-------------------!
        ! Loading variable: !
        !-------------------!

        ! get variable ID and number of dimension:
        ierr = nf90_inq_varid( fileid, varname, varid )
        call nf90check(ierr,'Error while getting identifiers of variable '//trim(varname))

        ! get variable fillvalue
        ierr = nf90_get_att( fileid, varid, fillvalname, fillval )
        call nf90check(ierr,'Warning: unable to get attribute "'//fillvalname//'" of variable '//trim(varname), kill=.false.)
        got_fillval = (ierr==NF90_NOERR)

        ! get variable shape:
        ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
        call nf90check(ierr,'Error while getting dimensions identifiers of variable '//trim(varname))
        do i = 1,3
            ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
            call nf90check(ierr,'Error while getting dimensions lengths and names of variable '//trim(varname))
        end do

        ! print report:
        print *, 'netCDF variable: '//trim(varname)
        print *, '  - Loaded dimension:'
        print *, '      lon:   ',dimname(1)
        print *, '      lat:   ',dimname(2)
        print *, '      litho: ',dimname(3)
        print *

        ! variable loading:
        if (dimlen(1)/=nlon .or. dimlen(2)/=nlat .or. dimlen(3)/=nlitho) then
            print *, 'Error while loading variable '//varname
            print *, 'Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),'x',dimlen(3)
            print *, 'while expected:          ',nlon,'x',nlat,'x',nlitho
            stop
        else
            do j=1,nlat
            do k=1,nlitho
                ierr = nf90_get_var( fileid, varid, litho_frac( k , 1+nlon*(j-1) : nlon*j ), start=(/1,j,k/), count=(/nlon,1,1/) )
                !         dimension unravelling --------------------^^^^^^^^^^^^^^^^^^^^^
                call nf90check(ierr,'Error while getting variable '//trim(varname))
            end do
            end do
        end if


        !----------------------------------------!
        ! Check consistency with reference axis: !
        !----------------------------------------!

        if (.not. (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0))) then ! if reference axis are defined
            ! Assume: 1st dim = x axis, 2nd dim = y axis
            ierr = nf90_inq_varid(fileid, dimname(1), varid)
            call nf90check(ierr,'Warning: unable to get identifier of variable '//trim(dimname(1)), kill=.false.)
            if (ierr==NF90_NOERR) then
                ierr = nf90_get_var(fileid, varid, loc_x_axis)
                call nf90check(ierr,'Error while getting variable '//trim(dimname(1)))
                ierr = nf90_inq_varid(fileid, dimname(2), varid)
                call nf90check(ierr,'Warning: unable to get identifier of variable '//trim(dimname(2)), kill=.false.)
                if (ierr==NF90_NOERR) then
                    ierr = nf90_get_var(fileid, varid, loc_y_axis)
                    call nf90check(ierr,'Error while getting variable '//trim(dimname(2)))
                    call check_axis('lithology file', loc_x_axis, loc_y_axis, ref_x_axis, ref_y_axis, ERROR_HANDLING_OPTION(1))
                else
                    print *, 'Cannot check axis consistency'
                end if
            else
                print *, 'Cannot check axis consistency'
            end if
        end if


        !---------------!
        ! File closing: !
        !---------------!

        ierr = nf90_close(fileid)
        call nf90check(ierr,'Error while closing input file '//trim(filename))


        !---------!
        ! Checks: !
        !---------!

        ! Missingpoints
        if (got_fillval) then
            call check_landcells('lithology fraction', fillval, areaclimber, ERROR_HANDLING_OPTION(2), var2D=litho_frac, axis=1)
        end if

        ! Fraction consistency
        call check_invalid('lithology fraction', areaclimber, ERROR_HANDLING_OPTION(4), var2D=litho_frac, axis=1)


    end if


    ! Close scratch file
    close(ID)


end subroutine
