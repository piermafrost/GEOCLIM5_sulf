module utils
! Miscellaneous file and string handling functions

implicit none
contains


    subroutine add_path(filename)
        character(len=*), intent(inout):: filename
        include 'path.inc' ! => 'geoclim_path' variable
        !
        ! test if absolute path given (file name starts with a '/')
        ! if not, add current GEOCLIM path to the name:
        if (filename(1:1) /= '/') then
            filename = geoclim_path//filename
        end if
    end subroutine


    !----------------------------------------------------------------------!


    subroutine read_comment(funit, ierr)
      ! read file until it finds uncommented line
      ! Commented lines are blank (empty) lines or lines whose first non-blank character is '#'
      integer, intent(in):: funit
      integer, intent(out), optional:: ierr
      character(len=1):: line

      line = '#'

      if (present(ierr)) then

        ! read until finding uncommented line or error is raised
        ierr = 0
        do while (line=='#' .and. ierr==0)
          read(unit=funit, fmt=*, iostat=ierr) line
        end do

        ! if no error raised, get back to the previous record (in that case, the previous line)
        if (ierr==0) then
          backspace(unit=funit)
        end if

      else

        ! read until finding uncommented line
        do while (line=='#')
          read(unit=funit, fmt=*) line
        end do

        ! get back to the previous record (in that case, the previous line)
        backspace(unit=funit)

      end if

    end subroutine


    !----------------------------------------------------------------------!


    subroutine set_error_handling_option()
      integer, dimension(5):: ERROR_HANDLING_OPTION
      common /error/ ERROR_HANDLING_OPTION
      character(len=30):: argchar
      integer:: arg, ierr

      ! option = -1: ask the user interactively (default)
      ! option = 0:  abort the program
      ! option = 1:  remove erratic points (set area=0). Not possible for axis or units error.
      ! option = 2:  ignore the issue and continue the execution
      ! option = 3:  replace invalid value. Only valid for negative runoff (> set 0) and null slope (> set minval) errors

      ! 1 Axis mismatch error
      ! ---------------------
      call get_command_argument(1, argchar, status=ierr)
      read(argchar, fmt=*, iostat=ierr) arg
      if (ierr==0 .and. (arg==0 .or. arg==2)) then
        ERROR_HANDLING_OPTION(1) = arg
      else
        ERROR_HANDLING_OPTION(1) = -1 ! default
      end if

      ! 2 Missing-value on continental cells (area>0)
      ! ---------------------------------------------
      call get_command_argument(2, argchar, status=ierr)
      read(argchar, fmt=*, iostat=ierr) arg
      if (ierr==0 .and. arg>=-1 .and. arg<=2) then
        ERROR_HANDLING_OPTION(2) = arg
      else
        ERROR_HANDLING_OPTION(2) = -1 ! default
      end if

      ! 3 Invalid value (negative runoff, null slope, ...)
      ! --------------------------------------------------
      ! read potential argument
      call get_command_argument(3, argchar, status=ierr)
      read(argchar, fmt=*, iostat=ierr) arg
      if (ierr==0 .and. arg>=-1 .and. arg<=3) then
        ERROR_HANDLING_OPTION(3) = arg
      else
        ERROR_HANDLING_OPTION(3) = -1 ! default
      end if

      ! 4 Sum of all lithological classes' fraction is not 1
      ! -----------------------------------------------------
      ! read potential argument
      call get_command_argument(4, argchar, status=ierr)
      read(argchar, fmt=*, iostat=ierr) arg
      if (ierr==0 .and. arg>=-1 .and. arg<=2) then
        ERROR_HANDLING_OPTION(4) = arg
      else
        ERROR_HANDLING_OPTION(4) = -1 ! default
      end if

      ! 5 Units not recognized error
      ! ----------------------------
      call get_command_argument(5, argchar, status=ierr)
      read(argchar, fmt=*, iostat=ierr) arg
      if (ierr==0 .and. (arg==0 .or. arg==2)) then
        ERROR_HANDLING_OPTION(5) = arg
      else
        ERROR_HANDLING_OPTION(5) = -1 ! default
      end if

    end subroutine
        

end module
