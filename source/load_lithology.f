subroutine load_lithology(ID)
    use netcdf
    use netcdf_io_module, only: nf90_check
    use io_module, only: check_landcells, check_axis, check_invalid, UNDEFINED_VALUE_CHAR, check_namelist_def
    use utils, only: add_path
    implicit none

    include 'combine_foam.inc'
    integer, intent(in):: ID
    character(len=500):: file_name
    character(len=100):: var_name, fillval_name
    integer:: fileid, varid, ierr
    integer, dimension(3):: dimids, dimlen
    character(len=50), dimension(3):: dimname
    double precision:: fillval
    logical:: got_fillval
    double precision, dimension(nlon):: loc_x_axis
    double precision, dimension(nlat):: loc_y_axis
    double precision, dimension(nlitho):: singlepixel_lithofrac ! unused variable, but present in namelist

    ! Namelist declaration
    namelist /LITHO_INFO/ file_name, var_name, fillval_name, singlepixel_lithofrac



    !----------------------------------------------!
    ! Get input information from main config file: !
    !----------------------------------------------!

    ! variables default value
    file_name    = UNDEFINED_VALUE_CHAR
    var_name     = UNDEFINED_VALUE_CHAR
    fillval_name = UNDEFINED_VALUE_CHAR

    rewind(unit=ID)
    ! <><><><><><><><><><><><><> !
    read(unit=ID, nml=LITHO_INFO)
    ! <><><><><><><><><><><><><> !

    call check_namelist_def('Error - in open_input_files.f: variable "file_name" from namelist "LITHO_INFO" was not given in'// &
                            ' config/IO_CONDITIONS', char_var=file_name)
    call check_namelist_def('Error - in open_input_files.f: variable "var_name" from namelist "LITHO_INFO" was not given in'// &
                            ' config/IO_CONDITIONS', char_var=var_name)
    call check_namelist_def('Error - in open_input_files.f: variable "fillval_name" from namelist "LITHO_INFO" was not given in'// &
                            ' config/IO_CONDITIONS', char_var=fillval_name)

    call add_path(file_name)


    ! ========== !


    print *
    print *
    print *, 'Load lithology map'
    print *, '------------------'
    print *


    !----------------!
    ! File openning: !
    !----------------!

    ierr = nf90_open(file_name, NF90_NOWRITE, fileid)
    call nf90_check(ierr, 'Error while openning input file '//file_name)


    !-------------------!
    ! Loading variable: !
    !-------------------!

    ! get variable ID and number of dimension:
    ierr = nf90_inq_varid( fileid, var_name, varid )
    call nf90_check(ierr, 'Error while getting identifiers of variable '//var_name)

    ! get variable fillvalue
    ierr = nf90_get_att( fileid, varid, fillval_name, fillval )
    call nf90_check(ierr, 'Warning: unable to get attribute "'//fillval_name//'" of variable '//trim(var_name), kill=.false.)
    got_fillval = (ierr==NF90_NOERR)

    ! get variable shape:
    ierr = nf90_inquire_variable( fileid, varid, dimids=dimids )
    call nf90_check(ierr, 'Error while getting dimensions identifiers of variable '//var_name)
    do i = 1,3
        ierr = nf90_inquire_dimension( fileid, dimids(i), dimname(i), dimlen(i) )
        call nf90_check(ierr, 'Error while getting dimensions lengths and names of variable '//var_name)
    end do

    ! print report:
    print *, 'netCDF variable: '//trim(var_name)
    print *, '  - Loaded dimension:'
    print *, '      lon:   ',dimname(1)
    print *, '      lat:   ',dimname(2)
    print *, '      litho: ',dimname(3)
    print *

    ! variable loading:
    if (dimlen(1)/=nlon .or. dimlen(2)/=nlat .or. dimlen(3)/=nlitho) then
        print *, 'Error while loading variable '//var_name
        print *, 'Inconsistent shape, get: ',dimlen(1),'x',dimlen(2),'x',dimlen(3)
        print *, 'while expected:          ',nlon,'x',nlat,'x',nlitho
        stop
    else
        do j=1,nlat
        do k=1,nlitho
            ierr = nf90_get_var( fileid, varid, litho_frac( k , 1+nlon*(j-1) : nlon*j ), start=(/1,j,k/), count=(/nlon,1,1/) )
            !         dimension unravelling --------------------^^^^^^^^^^^^^^^^^^^^^
            call nf90_check(ierr, 'Error while getting variable '//var_name)
        end do
        end do
    end if


    !----------------------------------------!
    ! Check consistency with reference axis: !
    !----------------------------------------!

    if (.not. (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0))) then ! if reference axis are defined
        ! Assume: 1st dim = x axis, 2nd dim = y axis
        ierr = nf90_inq_varid(fileid, dimname(1), varid)
        call nf90_check(ierr, 'Warning: unable to get identifier of variable '//dimname(1), kill=.false.)
        if (ierr==NF90_NOERR) then
            ierr = nf90_get_var(fileid, varid, loc_x_axis)
            call nf90_check(ierr, 'Error while getting variable '//dimname(1))
            ierr = nf90_inq_varid(fileid, dimname(2), varid)
            call nf90_check(ierr, 'Warning: unable to get identifier of variable '//dimname(2), kill=.false.)
            if (ierr==NF90_NOERR) then
                ierr = nf90_get_var(fileid, varid, loc_y_axis)
                call nf90_check(ierr, 'Error while getting variable '//dimname(2))
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
    call nf90_check(ierr, 'Error while closing input file '//file_name)


    !---------!
    ! Checks: !
    !---------!

    ! Missingpoints
    if (got_fillval) then
        call check_landcells('lithology fraction', fillval, areaclimber, ERROR_HANDLING_OPTION(2), var2D=litho_frac, axis=1)
    end if

    ! Fraction consistency
    call check_invalid('lithology fraction', areaclimber, ERROR_HANDLING_OPTION(4), var2D=litho_frac, axis=1)


end subroutine
