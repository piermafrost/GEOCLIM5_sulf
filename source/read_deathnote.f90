subroutine read_deathnote(t)
implicit none
include 'combine_foam.inc'
integer:: ierr

open(unit=666, status='old', action='read', file=killing_file_name, iostat=ierr)
!                           do not raise error if file do not exist ^^^^^^^^^^^

do while(ierr==0)
  check_run_name = ' '
  read(unit=666,fmt=*,iostat=ierr) check_run_name
  if (check_run_name==run_name) then
    tend = t - 10d6 ! to end iteration loop
    ageYprint = t - 10d6 ! to create restart
    icount_geogprint = -1 !
    icount_DS_pri = -1    ! don't write any output
  end if
end do

close(unit=666)

end subroutine
