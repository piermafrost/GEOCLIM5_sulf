subroutine open_ascii_files(ID)
    use utils, only: add_path, read_comment
    implicit none
    integer, intent(in):: ID
    include 'combine_foam.inc'
    character(len=1000):: line
    integer:: ierr
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                INPUT AND INITIALIZATION FILES               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !GEOCLIM_VAR_INIT_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (2  ,status='old',file=filename)
    
    !SPECIES_INIT_FILE:
    call read_comment(ID)
    if (coupling_ecogeo) then
        print *, 'ECOGEO module not available'
        stop
    else
        read(ID,*)
    end if
    
    !BIODIV_VAR_INIT_FILE:
    call read_comment(ID)
    if (coupling_ecogeo) then
        print *, 'ECOGEO module not available'
        stop
    else
        read(ID,*)
    end if
    
    
    !--------!
    
    
    !CONTINENTAL_INPUT_MODE
    call read_comment(ID)
    read(unit=ID, fmt=*) dummychar, land_input_mode
    
    if (land_input_mode=='ascii') then
    
        !GRID_AREA_FILE:
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        call add_path(filename)
        open (301,status='old',file=filename)
    
        !CONTINENTAL_AREA_FILE:
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        call add_path(filename)
        open (7  ,status='old',file=filename)
    
    
        !ccccccc +++++++++++++++++++++++++++++++++ 
        !       Climate under several CO2 levels
    
        !TEMPERATURE_FILE:
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        call add_path(filename)
        open (30 ,status='old',file=filename)
    
        !RUNOFF_FILE:
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        call add_path(filename)
        open (31 ,status='old',file=filename)
    
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        !call add_path(filename)
        !open (302 ,status='old',file=filename)
    
        call read_comment(ID)
        read(unit=ID,fmt=*) dummychar, filename
        !call add_path(filename)
        !open (303 ,status='old',file=filename)
    
        !GCM_INPUT_CONDITIONS_FILE (ignore)
        call read_comment(ID)
        read(unit=ID, fmt=*)
    
    
    elseif (land_input_mode == 'GCM') then
    
        ! Skip 6 ascii-mode entries
        do k = 1,6
            call read_comment(ID)
            read(unit=ID, fmt=*)
        end do
    
        !GCM_INPUT_CONDITIONS_FILE
        call read_comment(ID)
        read(unit=ID, fmt=*) dummychar, filename
        call add_path(filename)
        open (333, status='old', file=filename)
    
    
    else
    
        print *, 'ERROR: continental input mode (defined config/IO_CONDITIONS) must be "ascii" or "GCM"'
        stop
    
    end if
    
    
    
    !LITHOLOGICAL FRACTION (NETCDF FILE): 
    call read_comment(ID)
    ! => store info in scratch file (that will be read by 'load_lithology' subroutine)
    read(unit=ID,fmt='(A1000)') line
    open (304 ,status='scratch')
    write(304, fmt='(A1000)') line
    
    !VEGETATION_FILE:
    call read_comment(ID)
    if (coupling_veget==1) then
        read(unit=ID,fmt=*) dummychar, filename
        call add_path(filename)
        open (48 ,status='old',file=filename)
    else
        read(ID,*)
    end if
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                     COMBINE INPUT FILES                     !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !OCEANIC_TEMPERATURE_FILE
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    ! Check if got "parametric" instead of file name
    if (filename=='parametric') then
        ! if so, store that info in a scratch file
        open(unit=32, status='scratch', action='readwrite')
        write(unit=32, fmt=*) 'parametric'
        rewind(unit=32)
    else
        ! if not, open expected input ascii file
        call add_path(filename)
        open (32 ,status='old',file=filename)
    end if
    
    ! **OBSOLETE** (former HYPSO_FILE):
    call read_comment(ID)
    read(unit=ID,fmt=*)! dummychar, filename
    !call add_path(filename)
    !open (4  ,status='old',file=filename)
    
    
    !ccccccc +++++++++++++++++++++++++++++++++ 
    !       Earth physical dimensions
    
    !OCEANIC_BOX_VOLUME_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (33 ,status='old',file=filename)
    
    !OCEAN_SURFACE_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (34 ,status='old',file=filename)
    
    !SURF_SEDIM_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (37 ,status='old',file=filename)
    
    
    !ccccccc+++++++++++++++++++++++++++++++++
    !       definition of boxes
    
    !DEEP_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (35 ,status='old',file=filename)
    
    !SEDIM_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (36 ,status='old',file=filename)
    
    !BOX_TEMPERATURE_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (38 ,status='old',file=filename)
    
    !APP_CONT_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (39 ,status='old',file=filename)
    
    !THERMOCLINE_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (40 ,status='old',file=filename)
    
    !SURFACE_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (41 ,status='old',file=filename)
    
    !EXCHANGE_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (42 ,status='old',file=filename)
    
    !
    ! unit=43: obsolete file: INPUT/indice_part.dat
    !
    
    !FSINK_INORG_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (44 ,status='old',file=filename)
    
    !FSINK_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (45 ,status='old',file=filename)
    
    !EPICONT_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (46 ,status='old',file=filename)
    
    !POLAR_BOX_FILE:
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (47 ,status='old',file=filename)
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                  RESTART AND OFFLINE FILES                  !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !ccccccc+++++++++++++++++++++++++++++++++ 
    !       geoclim restart files:
    
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    open (10 ,status='REPLACE',file=trim(output_path)//trim(filename)//run_name)
    
    
    
    if (coupling_ecogeo) then
        print *, 'ECOGEO module not available'
        stop
    else
    
        do i = 1,12
            call read_comment(ID)
            read(ID,*)
        end do
    
    end if
    
    
    
    
end subroutine
