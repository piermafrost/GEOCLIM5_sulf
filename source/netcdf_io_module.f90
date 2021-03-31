module netcdf_io_module
implicit none

contains


 subroutine create( filename, ID , mode )
 use netcdf

  character(len=*), intent(in):: filename
  integer, intent(out):: ID
  integer, intent(in), optional:: mode
  integer:: ierr

  if (present(mode)) then
    ierr = nf90_create( filename, mode, ID )
  else
    ierr = nf90_create( filename, NF90_CLOBBER, ID )
  end if
  call isitok(ierr,'Error while creating output file '//filename)

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine open_file( filename, ID , mode )
 use netcdf

  character(len=*), intent(in):: filename
  integer, intent(out):: ID
  integer, intent(in), optional:: mode
  integer:: ierr

  if (present(mode)) then
    ierr = nf90_open( filename, mode, ID )
  else
    ierr = nf90_open( filename, NF90_NOWRITE, ID )
  end if
  call isitok(ierr,'Error while openning file '//filename)

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine enddef( ID )
 use netcdf

  integer, intent(in):: ID
  character(len=30):: num
  integer:: ierr

  ierr = nf90_enddef( ID )
  write(num,fmt="(I0)") ID
  call isitok(ierr,'Error while end of definition of output file #'//trim(num))

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine redef( ID )
 use netcdf

  integer, intent(in):: ID
  character(len=30):: num
  integer:: ierr

  ierr = nf90_redef( ID )
  write(num,fmt="(I0)") ID
  call isitok(ierr,'Error while returning into define mode of output file #'//trim(num))

 end subroutine


!-------------------------------------------------------------------------------------------------!

 subroutine close_file( ID )
 use netcdf

  integer, intent(in):: ID
  character(len=30):: num
  integer:: ierr

  ierr = nf90_close( ID )
  write(num,fmt="(I0)") ID
  call isitok(ierr,'Error while closing output file #'//trim(num))

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine def_dim( ID, dimname , length, dimid )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), dimension(:), intent(in):: dimname
  integer, dimension(:), intent(in):: length
  integer, dimension(:), intent(out):: dimid
  integer:: ierr, N, i

  N = size(dimname,1)

  !===============================================================!
  do i=1,N
    ierr = nf90_def_dim( ID, dimname(i) , length(i), dimid(i) )
    call isitok(ierr,'Error while defining dimension '//dimname(i))
  end do
  !===============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine def_var( ID, varname, vartype, dimid, varid )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), dimension(:), intent(in):: varname
  integer, dimension(:), intent(in):: vartype, dimid
  integer, dimension(:), intent(out):: varid
  integer, allocatable:: vec_types(:)
  integer:: ierr, N, i

  N = size(varid,1)
  allocate(vec_types(N))

  if (size(vartype,1)==1) then
    vec_types = vartype(1)
  else
    vec_types = vartype
  end if


  !==================================================================!
  do i=1,N
    ierr=nf90_def_var( ID, varname(i), vec_types(i), dimid, varid(i) )
    call isitok(ierr,'Error while defining variable '//varname(i))
  end do
  !==================================================================!


 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_att_text( ID, varid, attname, attribute )
 use netcdf

  integer, intent(in):: ID ! file identifier
  integer, dimension(:), intent(in):: varid
  character(len=*), dimension(:), intent(in):: attname, attribute
!  character(len=:), dimension(:), allocatable:: vec_attname, vec_attribute
  character(len=200), dimension(:), allocatable:: vec_attname, vec_attribute
  integer:: ierr, l, l2, N, i
  character(len=30):: num

  N = size(varid,1)
!  l = len(attname)
!  l2 = len(attribute)
!  allocate( character(len=l)::  vec_attname(N)   )
!  allocate( character(len=l2):: vec_attribute(N) )
! This syntax is not supported by ifort10.1.021, character variables length must be imposed
  allocate( vec_attname(N)   )
  allocate( vec_attribute(N) )

  if (size(attname,1)==1) then
    vec_attname = attname(1)
  else
    vec_attname = attname
  end if

  if (size(attribute,1)==1) then
    vec_attribute = attribute(1)
  else
    vec_attribute = attribute
  end if

  !======================================================================================!
  do i=1,N
    ierr=nf90_put_att( ID, varid(i) , vec_attname(i), vec_attribute(i) )
    write(num,fmt="(I0)") i
    call isitok(ierr,'Error while putting variable #'//trim(num)//' attribute: '//vec_attname(i))
  end do
  !======================================================================================!


 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_att_real( ID, varid, attname, attribute )
 use netcdf

  integer, intent(in):: ID ! file identifier
  integer, dimension(:), intent(in):: varid
  character(len=*), dimension(:), intent(in):: attname
!  character(len=:), dimension(:), allocatable:: vec_attname
  character(len=200), dimension(:), allocatable:: vec_attname
  real, dimension(:), intent(in):: attribute
  real, dimension(:), allocatable:: vec_attribute
  integer:: ierr, N, l, i
  character(len=30):: num

  N = size(varid,1)
  l = len(attname)
!  allocate( character(len=l):: vec_attname(N) )
! This syntax is not supported by ifort10.1.021, character variables length must be imposed
  allocate( vec_attname(N) )
  allocate( vec_attribute(N) )

  if (size(attname,1)==1) then
    vec_attname = attname(1)
  else
    vec_attname = attname
  end if

  if (size(attribute,1)==1) then
    vec_attribute = attribute(1)
  else
    vec_attribute = attribute
  end if

  !======================================================================================!
  do i=1,N
    ierr=nf90_put_att( ID, varid(i) , vec_attname(i), vec_attribute(i) )
    write(num,fmt="(I0)") i
    call isitok(ierr,'Error while putting variable #'//trim(num)//' attribute: '//vec_attname(i))
  end do
  !======================================================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_att_dble( ID, varid, attname, attribute )
 use netcdf

  integer, intent(in):: ID ! file identifier
  integer, dimension(:), intent(in):: varid
  character(len=*), dimension(:), intent(in):: attname
!  character(len=:), dimension(:), allocatable:: vec_attname
  character(len=200), dimension(:), allocatable:: vec_attname
  double precision, dimension(:), intent(in):: attribute
  double precision, dimension(:), allocatable:: vec_attribute
  integer:: ierr, N, l, i
  character(len=30):: num

  N = size(varid,1)
  l = len(attname)
!  allocate( character(len=l):: vec_attname(N) )
! This syntax is not supported by ifort10.1.021, character variables length must be imposed
  allocate( vec_attname(N) )
  allocate( vec_attribute(N) )

  if (size(attname,1)==1) then
    vec_attname = attname(1)
  else
    vec_attname = attname
  end if

  if (size(attribute,1)==1) then
    vec_attribute = attribute(1)
  else
    vec_attribute = attribute
  end if

  !======================================================================================!
  do i=1,N
    ierr=nf90_put_att( ID, varid(i) , vec_attname(i), vec_attribute(i) )
    write(num,fmt="(I0)") i
    call isitok(ierr,'Error while putting variable #'//trim(num)//' attribute: '//vec_attname(i))
  end do
  !======================================================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_att_int( ID, varid, attname, attribute )
 use netcdf

  integer, intent(in):: ID ! file identifier
  integer, dimension(:), intent(in):: varid
  character(len=*), dimension(:), intent(in):: attname
!  character(len=:), dimension(:), allocatable:: vec_attname
  character(len=200), dimension(:), allocatable:: vec_attname
  integer, dimension(:), intent(in):: attribute
  integer, dimension(:), allocatable:: vec_attribute
  integer:: ierr, N, l, i
  character(len=30):: num

  N = size(varid,1)
  l = len(attname)
!  allocate( character(len=l):: vec_attname(N) )
! This syntax is not supported by ifort10.1.021, character variables length must be imposed
  allocate( vec_attname(N) )
  allocate( vec_attribute(N) )

  if (size(attname,1)==1) then
    vec_attname = attname(1)
  else
    vec_attname = attname
  end if

  if (size(attribute,1)==1) then
    vec_attribute = attribute(1)
  else
    vec_attribute = attribute
  end if

  !======================================================================================!
  do i=1,N
    ierr=nf90_put_att( ID, varid(i) , vec_attname(i), vec_attribute(i) )
    write(num,fmt="(I0)") i
    call isitok(ierr,'Error while putting variable #'//trim(num)//' attribute: '//vec_attname(i))
  end do
  !======================================================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_real1D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  real, dimension(:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  real, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_real2D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  real, dimension(:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  real, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_real3D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  real, dimension(:,:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  real, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_real4D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  real, dimension(:,:,:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  real, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_dble1D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  double precision, dimension(:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  double precision, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_dble2D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  double precision, dimension(:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  double precision, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_dble3D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  double precision, dimension(:,:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  double precision, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine get_var_dble4D( ID, varname, VAR, units, missval_name ,missval, stt, cnt )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), intent(in):: varname
  double precision, dimension(:,:,:,:), intent(out):: VAR
  character(len=*), intent(out), optional:: units
  character(len=*), intent(in), optional:: missval_name
  double precision, intent(out), optional:: missval
  integer, intent(in), dimension(:), optional:: stt, cnt

  character(len=50):: missval_n
  integer:: varid, ierr, i, l, n

  ierr = nf90_inq_varid( ID, varname, varid )
  call isitok(ierr,'Error while getting '//varname//' variable identifiers')
  if (present(stt) .and. present(cnt)) then
    ierr = nf90_get_var( ID, varid, VAR , start=stt, count=cnt)
  else
    ierr = nf90_get_var( ID, varid, VAR )
  end if
  call isitok(ierr,'Error while getting variable '//varname)

  if (present(units)) then
    ierr = nf90_get_att( ID, varid, 'units', units)
    call isitok(ierr,'Error while getting variable '//varname//' attribute "units"')
  end if

  if (present(missval)) then
    if (present(missval_name)) then
      l=len(missval_name)
      missval_n(1:l) = missval_name
    else
      l=10
      missval_n(1:l) = '_Fillvalue'
    end if
    ierr = nf90_get_att( ID, varid, missval_n(1:l), missval)
    call isitok(ierr,'Error while getting variable '//varname//' attribute '//missval_n(1:l))
  end if

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_real1D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  real, dimension(:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_real2D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  real, dimension(:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_real3D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  real, dimension(:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_real4D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  real, dimension(:,:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_dble1D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  double precision, dimension(:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_dble2D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  double precision, dimension(:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_dble3D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  double precision, dimension(:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_dble4D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  double precision, dimension(:,:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_int1D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  integer, dimension(:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_int2D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  integer, dimension(:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_int3D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  integer, dimension(:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine put_var_int4D( ID, varid, VAR, begin, length )
 use netcdf

  integer, intent(in):: ID, varid
  integer, dimension(:,:,:,:), intent(in):: VAR
  integer, dimension(:), optional:: begin, length
  integer:: ierr
  character(len=30):: num

  !==============================================================!
  if ( present(begin) .and. present(length) ) then
    ierr=nf90_put_var( ID, varid ,VAR, start=begin, count=length )
  else
    ierr=nf90_put_var( ID, varid ,VAR )
  end if
  write(num,fmt="(I0)") varid
  call isitok(ierr,'Error while putting variable #'//trim(num))
  !==============================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine inquire_dim( ID, dimname, dimid )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), dimension(:), intent(in):: dimname
  integer, dimension(:), intent(out):: dimid
  integer:: ierr, N, i

  N = size(dimname,1)

  !===============================================================================!
  do i=1,N
    ierr = nf90_inq_dimid( ID, dimname(i) , dimid(i) )
    call isitok(ierr,'Error while inquiring dimension '//dimname(i)//' identifier')
  end do
  !===============================================================================!

 end subroutine


!-------------------------------------------------------------------------------------------------!


 subroutine inquire_var( ID, varname, varid )
 use netcdf

  integer, intent(in):: ID ! file identifier
  character(len=*), dimension(:), intent(in):: varname
  integer, dimension(:), intent(out):: varid
  integer:: ierr, N, i

  N = size(varid,1)

  !==============================================================================!
  do i=1,N
    ierr=nf90_inq_varid( ID, varname(i), varid(i) )
    call isitok(ierr,'Error while inquiring variable '//varname(i)//' identifier')
  end do
  !==============================================================================!


 end subroutine


!-------------------------------------------------------------------------------------------------!


  subroutine check_fillvalue( ID, varid, fillvalname, fillval )
  use netcdf

   integer, intent(in):: ID ! file identifier
   integer, dimension(:), intent(in):: varid
   character(len=*), dimension(:), intent(in):: fillvalname
   real, dimension(:), intent(inout):: fillval
   real:: fillvalue
   integer:: ierr, N, i
   character(len=30):: num

   N = size(varid,1)

   !=======================================================================================!
   do i=1,N
     ierr=nf90_get_att( ID, varid(i), fillvalname(i), fillvalue )
     write(num,fmt="(I0)") varid(i)
     call isitok(ierr,'Error while checking variable #'//trim(num)//' attribute '//fillvalname(i))
       if ( ierr/=NF90_NOERR ) then
         print *, 'CANNOT TEST MISSING-VALUE CONSISTENCY. PROGRAM STOPPED'
         stop
       else
         if ( fillvalue /= fillval(i) ) then
           print *, 'WARNING: variable #'//trim(num)//' missing-value is not consistent with'
           print *, 'existing file. It will be changed to keep existing file value' 
           fillval(i) = fillvalue
         end if
       end if
   end do
   !======================================================================================!


  end subroutine


!-------------------------------------------------------------------------------------------------!

 subroutine isitok(ierr,message)
 use netcdf
   integer, intent(in):: ierr
   character(len=*), intent(in):: message
   if (ierr/=NF90_NOERR) then
     print *
     print *, message
     print *, nf90_strerror(ierr)
     stop
   end if
 end subroutine




end module
