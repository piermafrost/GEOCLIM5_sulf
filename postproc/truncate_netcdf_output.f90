program truncate
use netcdf
implicit none

integer, parameter:: n_dim_max=4
character(len=100), dimension(2), parameter:: timenamelist = (/'time','time_counter'/)
integer, dimension(n_dim_max):: starts, length
character(len=100):: fname, dimname(n_dim_max), varname, attname, dummychar
character(len=2):: num
real:: timestart, timestop
real, allocatable, dimension(:):: t
real, allocatable, dimension(:,:,:,:):: VAR
integer:: ifileid, ofileid, ivarid, dimvarid(n_dim_max), Odimids(n_dim_max), Vdimids(n_dim_max), ierr
integer, allocatable:: Ovarids(:)
integer:: i0, Ndim, Nvar, Natt, n, i, j, l, TIMEDIM, stt(n_dim_max), ctt(n_dim_max), sttO(n_dim_max)
logical:: loop
integer:: nb_arg, iargc
external:: iargc


nb_arg = iargc()
if ( nb_arg > 3) then
   print *, "WARNING : too many input arguments (three at most expected), only the first three will be considered"
end if
!
if (nb_arg>=1) then
  call getarg(1,fname)
else
  write(*,fmt='(A)',advance='no') 'Input file name: '
  read *, fname
end if
!
if (nb_arg>=2) then
  call getarg(2,dummychar)
  read(dummychar,fmt=*) timestart
else
  write(*,fmt='(A)',advance='no') 'Truncating starting time: '
  read *, timestart
end if
!
if (nb_arg>=3) then
  call getarg(3,dummychar)
  read(dummychar,fmt=*) timestop
else
  write(*,fmt='(A)',advance='no') 'Truncating ending time: '
  read *, timestop
end if

!-----------------------------------------------------!

ierr = nf90_open(fname, NF90_NOWRITE, ifileid)
call nf90check(ierr,'Error while openning file '//fname)

ierr = nf90_inquire( ifileid, nDimensions=Ndim )
call nf90check(ierr,'Error while getting input file number of dimension')

if (Ndim > n_dim_max) then
  write(num,fmt='(I2)') n_dim_max
  print *, 'ERROR, too many dimensions. Maximum supported number of dimensions is '//num
  call abort
end if

starts = 1
length = 1
TIMEDIM=0
do i=1,Ndim
  ierr = nf90_inquire_dimension( ifileid, i, name=dimname(i), len=l )
  write(num,fmt='(I2)') i
  call nf90check(ierr,'Error while getting dimension #'//num//' name')
  length(i) = l
  dummychar = dimname(i)
  if (any( dummychar == timenamelist )) then
    TIMEDIM=i
  end if
end do

if (TIMEDIM==0) then
  print *, 'ERROR: no time dimension found'
  call abort
end if

do i=1,Ndim
  ierr = nf90_inq_varid(ifileid,dimname(i),dimvarid(i))
  write(num,fmt='(I2)') i
  call nf90check(ierr,'Error while getting dimension #'//num//' variable ID')
end do

allocate(t(length(TIMEDIM)))

ierr = nf90_get_var(ifileid,dimvarid(TIMEDIM),t)
call nf90check(ierr,'Error while getting time variable')

i0 = 1
i = 1
loop=.true.
do while (loop)
  i = i+1
  if ( i > length(TIMEDIM) ) then
    loop = .false.
    i = i-1
  else
    if ( t(i) > timestop ) then
      loop = .false.
      i = i-1
    else
      if ( t(i0) < timestart ) i0 = i0+1
    end if
  end if
end do
starts(TIMEDIM) = i0
length(TIMEDIM) = i-i0+1



allocate(VAR(length(1),length(2),length(3),length(4)))



ierr = nf90_inquire( ifileid, nVariables=Nvar )
call nf90check(ierr,'Error while getting input file number of Variables')

allocate(Ovarids(Nvar))


ierr = nf90_create( 'TRUNCATED_'//fname, NF90_CLOBBER, ofileid  )
call nf90check(ierr,'Error while creating output file TRUNCATED_'//fname)

 ! Dimensions:
do i=1,Ndim
  if (i==TIMEDIM) then
    l=NF90_UNLIMITED
  else
    l=length(i)
  end if
  ierr = nf90_def_dim( ofileid, dimname(i) , len=l  , dimid=Odimids(i) )
  call nf90check(ierr,'Error while defining dimension '//dimname(i))
end do


do i=1,Nvar

  ierr = nf90_inquire_variable( ifileid, varid=i,  ndims=n )
  write(num,fmt='(I2)') i
  call nf90check(ierr,'Error while inquiring variable #'//num//' number of dimensions')
  ierr = nf90_inquire_variable( ifileid, varid=i, name=varname, dimids=Vdimids(1:n), nAtts=Natt )
  call nf90check(ierr,'Error while inquiring variable #'//num)

  ierr = nf90_def_var( ofileid, varname, NF90_FLOAT, Vdimids(1:n), Ovarids(i) )
  call nf90check(ierr,'Error while defining variable'//varname)

  do j=1,Natt
    ierr = nf90_inq_attname( ifileid, i, attnum=j, name=attname)
    write(num,fmt='(I2)') j
    call nf90check(ierr,'Error while getting attribute #'//num//' name of variable '//varname)
    ierr = nf90_copy_att( ifileid, i, attname, ofileid, Ovarids(i) )
    call nf90check(ierr,'Error while copying attribute '//attname//' of variable '//varname)
  end do

end do

ierr = nf90_enddef(ofileid)

stt = 1
sttO = 1
ctt = 1

do i=1,Nvar

  ierr = nf90_inquire_variable( ifileid, varid=i,  ndims=n )
  write(num,fmt='(I2)') i
  call nf90check(ierr,'Error while inquiring variable #'//num//' number of dimensions')
  ierr = nf90_inquire_variable( ifileid, varid=i, name=varname, dimids=Vdimids(1:n), nAtts=Natt )
  call nf90check(ierr,'Error while inquiring variable #'//num)

  do j=1,n
    stt(j) = starts(Vdimids(j))
    ctt(j) = length(Vdimids(j))
  end do

  ierr = nf90_get_var( ifileid, varid=i, values=VAR, start=stt(1:n), count=ctt(1:n) ) 
  call nf90check(ierr,'Error while getting variable '//varname)
  ierr = nf90_put_var( ofileid, Ovarids(i), VAR, start=sttO(1:n), count=ctt(1:n) )
  call nf90check(ierr,'Error while putting variable '//varname)

end do

ierr = nf90_close(ifileid)
ierr = nf90_close(ofileid)

print *, 'File written: TRUNCATED_'//fname



end program truncate



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine nf90check(ierr,message)
use netcdf
  integer, intent(in):: ierr
  character(len=*), intent(in):: message
  if (ierr/=NF90_NOERR) then
    print *, message
    print *, nf90_strerror(ierr)
    call abort
  end if
end subroutine

