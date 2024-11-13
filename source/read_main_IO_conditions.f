subroutine read_main_IO_conditions(ID)
    use utils, only: add_path
    use io_module, only: UNDEFINED_VALUE_CHAR, check_namelist_def
    implicit none
    integer, intent(in):: ID
    include 'combine_foam.inc'

    ! Local variables
    character(len=500):: phys_cond_file

    ! Namelist declaration
    namelist /MAIN_INFO/ run_name, output_directory, phys_cond_file, killing_file

    ! Default values of namelist variables:
    run_name         = UNDEFINED_VALUE_CHAR  ! Name of run
    output_directory = UNDEFINED_VALUE_CHAR  ! Directory where to write the outputs
    phys_cond_file   = UNDEFINED_VALUE_CHAR  ! Configuration file with "physical" conditions
    killing_file     = UNDEFINED_VALUE_CHAR  ! File for killing signal (e.g., deathnote.txt), unit=666


    ! Read namelist variables in configuration file
    ! <><><><><><><><><><><><><><><><><><><><><><> !
    read(unit=ID, nml=MAIN_INFO)
    ! <><><><><><><><><><><><><><><><><><><><><><> !


    ! Check that all variables were set
    call check_namelist_def('Error - in read_main_IO_conditions.f: variable "run_name" was not given in config/IO_CONDITIONS', &
                            char_var=run_name)
    call check_namelist_def('Error - in read_main_IO_conditions.f: variable "output_directory" was not given in'// &
                            ' config/IO_CONDITIONS', char_var=output_directory)
    call check_namelist_def('Error - in read_main_IO_conditions.f: variable "phys_cond_file" was not given in'// &
                            ' config/IO_CONDITIONS', char_var=phys_cond_file)
    call check_namelist_def('Error - in read_main_IO_conditions.f: variable "killing_file" was not given in config/IO_CONDITIONS', &
                            char_var=killing_file)

    ! Add GEOCLIM root path to file names:
    call add_path(output_directory)
    call add_path(phys_cond_file)
    call add_path(killing_file)

    ! Open physical_conditions file (e.g., cond_p20.dat):
    open (3  , status='old', action='read', file=phys_cond_file)
    !    *******************************************************

end subroutine

