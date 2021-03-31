subroutine read_main_IO_conditions(ID)
    use utils, only: add_path, read_comment
    implicit none
    integer, intent(in):: ID
    include 'combine_foam.inc'

    ! Name of run
    call read_comment(ID)
    read(ID,*) dummychar, run_name

    ! Directory to write the outputs in
    call read_comment(ID)
    read(ID,*) dummychar, output_path
    call add_path(output_path)

    ! PHYSICAL_CONDITIONS_FILE (cond_p20):
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, filename
    call add_path(filename)
    open (3  ,status='old',file=filename)

    ! File for killing signal (deathnote.txt), unit=666
    call read_comment(ID)
    read(unit=ID,fmt=*) dummychar, killing_file_name
    call add_path(killing_file_name)

end
