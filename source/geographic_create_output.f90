module geographic_create_output_mod
implicit none

contains

subroutine geographic_create_output( ID,output_path,run_name,lon,lat,area,slope,litho_frac, fname,fnum,varname,missval )
  use io_module, only: read_io_condition
  use netcdf_io_module, only: create, def_dim, def_var, put_att_text, put_att_real, enddef, redef, put_var_real1D, &
                              put_var_int1D, put_var_real2D, put_var_real3D, close_file
  use netcdf
  include 'coupler.inc'
  integer, parameter:: nvar=15
  integer, intent(in):: ID
  character(len=*), intent(in):: output_path, run_name
  double precision, dimension(:), intent(in):: lon, lat, area
  double precision, dimension(:), intent(inout):: slope
  double precision, dimension(:,:), intent(inout):: litho_frac
  !
  integer, intent(out), dimension(nvar):: fnum
  character(len=*), intent(out), dimension(nvar) :: fname, varname
  double precision, intent(out), dimension(nvar):: missval
  character(len=100), dimension(nvar) :: varunits, missvalname
  character(len=500), dimension(nvar) :: var_longname
  !
  integer, dimension(nvar) :: fid, varid
  !
  integer, parameter:: maxlgth=500
  character(len=maxlgth):: dummy
  integer:: nbasin, nlon, nlat, nlit, dimid(nvar,4)
  integer:: i, j, k

  ! Default name and value of missing-values (if unspecified in text input files)
  character(len=*), parameter:: defmissvalname = '_FillValue'
  double precision, parameter:: defmissval = 9.96921e+36

  !===================================== OUTPUT VARIABLES LIST: ===================================================================!
  !  OUTPUT: lon, lat, litho time, area, litho_frac temperature, runoff, slope, sil_weathering, weighted_weathering                !
  !  var #:  1    2    3     4     5     6          7            8       9      10              11                                 !
  !          unweighted_weathering, ker_weathering, bio_Corg_export, phos_weathering                                               !
  !          12                     13              14               15                                                            !
  !================================================================================================================================!




  ! default missing-value
  missvalname = defmissvalname
  missval = defmissval

  ! get dimension size:
  nlon = size(lon)
  nlat = size(lat)
  nlit = size(litho_frac,1)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%                                   getting output file names:                                   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


  ! read variables attributes (#1 to #nvar)
  do i=1,4
    call read_io_condition(ID, fname(i) , varname(i) , varunits(i) )
  end do
  do i=5,nvar
    call read_io_condition(ID, fname(i) , varname(i) , varunits(i) , missvalname(i)  , missval(i) , var_longname(i) )
  end do


  do i=2,nvar
  ! replace "-" by previous line's file name
    dummy = fname(i)
    if (dummy(1:1) == '-') then
      dummy=''
      j = i
      do while (dummy=='')
        j = j-1
        dummy = fname(j)
      end do
      fname(i) = dummy
    else
      if (dummy(1:1) /= '#') then
        ! add output directory path and extension:
        fname(i) = trim(output_path)//trim(dummy)//trim(run_name)//'.nc'
      else
        fname(i) = ''
      end if
    end if
  end do



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%                                    output file definition:                                     %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


  ! Check the number of different files:
  !-------------------------------------

  fnum(6) = 6 ! the first 5 variables are dimensions variables + area, and don't need a file
  do k = 6,nvar  ! loop on all file names
    i = 1
    do while ( fname(k) /= fname(i) ) ! loop on all previous file names
      i = i + 1
    end do
    fnum(k) = i
  end do
  do k = 1,nvar
    if (fname(k)=='') then
      ! "#" file name => don't write this variable
      ! => put 0, so that it will not be created
      fnum(k) = 0
    end if
  end do


  ! Output file creation:
  !---------------------

  k = 0
  do i = 6,nvar
    if ( fnum(i) > k ) then ! if "new" file

      call create(fname(i), fid(i))
      k = fnum(i)

      ! GOBAL ATTRIBUTES
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'title'/),                       (/'Main continental ouputs'/))
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'run_name'/),                    (/trim(run_name)/))
      if (coupling_dynsoil) then
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'silicate_weathering_model'/), (/'DynSoil'/))
        if (use_dynsoil_steady_state) then
          call put_att_text(fid(i), (/NF90_GLOBAL/), (/'DynSoil_mode'/),            (/'steady-state'/))
        else
          call put_att_text(fid(i), (/NF90_GLOBAL/), (/'DynSoil_mode'/),            (/'dynamic'/))
        end if
      else
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'silicate_weathering_model'/), (/'old'/))
      endif
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'CO2_interpolation'/),           (/interpolation_mode/))

    elseif(fnum(i)/=0) then
      fid(i) = fid(fnum(i))

    else
      fid(i) = 0

    end if
  end do 


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Dimensions and dimension variables (definition, attribute, and putting): !
  ! [need to be defined as many time as different files]                     !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


  ! horizontal dimensions:
  !-----------------------

  k = 0
  do i = 6,nvar ! loop on all non-dimension variables nor area (as there are all at least 2D)
    if ( fnum(i) > k ) then
      call def_dim( fid(i), varname(1:2), (/nlon,nlat/) , dimid(i,1:2)                )
      call def_var( fid(i), varname(1:1), (/NF90_FLOAT/), dimid(i,1:1), varid(1:1)    )
      call def_var( fid(i), varname(2:2), (/NF90_FLOAT/), dimid(i,2:2), varid(2:2)    )
      ! dimension variables attributes:
      call put_att_text( fid(i), varid(1:1), (/'axis'/), (/'X'/)                      )
      call put_att_text( fid(i), varid(1:1), (/'nav_model'/), (/'Default grid'/)      )
      call put_att_text( fid(i), varid(2:2), (/'axis'/), (/'Y'/)                      )
      call put_att_text( fid(i), varid(2:2), (/'nav_model'/), (/'Default grid'/)      )
      ! Variable attributes:
      call put_att_text( fid(i), varid(1:2), (/'name'/) , varname(1:2)                )
      call put_att_text( fid(i), varid(1:2), (/'units'/), varunits(1:2)               )
      call enddef( fid(i) )
      ! Put variable:
      call put_var_real1D( fid(i), varid(1), real(lon)  )
      call put_var_real1D( fid(i), varid(2), real(lat)  )
      call redef( fid(i) )
      !
      k = fnum(i)
    elseif(fnum(i)/=0) then ! if old file
      dimid(i,1:2) = dimid(fnum(i),1:2)
    end if
  end do


  ! lithology dimension variable:
  !------------------------------

  k = 0
  i = 6 ! loop on all lithology-defined variables
  if ( fnum(i) > k ) then
    call def_dim( fid(i), varname(3:3), (/nlit/) , dimid(i,3:3)                     )
    call def_var( fid(i), varname(3:3), (/NF90_INT/), dimid(i,3:3), varid(3:3)      )
    ! Variable attributes:
    call put_att_text( fid(i), varid(3:3), (/'name'/) , varname(3:3)                )
    call put_att_text( fid(i), varid(3:3), (/'units'/), varunits(3:3)               )
    call enddef( fid(i) )
    ! Put variable:
    call put_var_int1D( fid(i), varid(3), (/(j,j=1,nlit)/) )
    call redef( fid(i) )
    !
    k = fnum(i)
  elseif(fnum(i)/=0) then ! if old file
    dimid(i,3) = dimid(fnum(i),3)
  end if
  do i = 11,12 ! loop on all lithology-defined variables
    if ( fnum(i) > k ) then
      call def_dim( fid(i), varname(3:3), (/nlit/) , dimid(i,3:3)                     )
      call def_var( fid(i), varname(3:3), (/NF90_INT/), dimid(i,3:3), varid(3:3)      )
      ! Variable attributes:
      call put_att_text( fid(i), varid(3:3), (/'name'/) , varname(3:3)                )
      call put_att_text( fid(i), varid(3:3), (/'units'/), varunits(3:3)               )
      call enddef( fid(i) )
      ! Put variable:
      call put_var_int1D( fid(i), varid(3), (/(j,j=1,nlit)/) )
      call redef( fid(i) )
      !
      k = fnum(i)
    elseif(fnum(i)/=0) then ! if old file
      dimid(i,3) = dimid(fnum(i),3)
    end if
  end do


  ! time dimensions and area variables:
  !------------------------------------

  k = 0
  do i = 7,nvar ! loop on all non-dimension variables nor area and litho_frac (as there are all time-defined)
    if ( fnum(i) > k ) then
      if (i/=9) then ! exclude slope variable
        call def_dim( fid(i), varname(4:4), (/NF90_UNLIMITED/) , dimid(i,4:4)           )
        call def_var( fid(i), varname(4:4), (/NF90_FLOAT/), dimid(i,4:4), varid(4:4)    )
        call def_var( fid(i), varname(5:5), (/NF90_FLOAT/), dimid(i,1:2), varid(5:5)    )
        ! dimension variables attributes:
        call put_att_text( fid(i), varid(4:4), (/'axis'/), (/'T'/)                      )
        ! Variable attributes:
        call put_att_text( fid(i), varid(4:5), (/'name'/) , varname(1:4)                )
        call put_att_text( fid(i), varid(4:5), (/'units'/), varunits(1:4)               )
        call put_att_real( fid(i), varid(5:5), missvalname(4:4) , real(missval(4:4))    )
        if (.not. (var_longname(5) == '')) then
          call put_att_text(fid(i), varid(5:5), (/'long_name'/), var_longname(5:5) ) ! area
        end if
        call enddef( fid(i) )
        ! Put variable:
        call put_var_real2D( fid(i), varid(5), real(1e12*reshape(area, shape=(/nlon,nlat/))) )
        call redef( fid(i) )
        !
        do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,4) = dimid(i,4) ; end do
        !
        k = fnum(i)
      end if
    elseif(fnum(i)/=0) then ! if old file
      dimid(i,4) = dimid(fnum(i),4)
    end if
  end do


  !++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Other variables [need to be defined one time only] !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++!


  ! Variables definition:
  !----------------------

  ! 3D lon+lat+litho variables (litho_frac)
  i = 6
  if (fnum(i)==0) then ! => not writing this variable
    fid(i) = 0
    varid(i) = 0
  else
    call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,1:3), varid(i:i) )
  end if

  ! 3D lon+lat+time variables
  do i = 7,8
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:2),dimid(i,4)/), varid(i:i) )
    end if
  end do

  ! 2D variables (slope)
  i = 9
  if (fnum(i)==0) then ! => not writing this variable
    fid(i) = 0
    varid(i) = 0
  else
    call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,1:2), varid(i:i) )
  end if

  ! 3D lon+lat+time variables (all sil wth)
  i = 10
  if (fnum(i)==0) then ! => not writing this variable
    fid(i) = 0
    varid(i) = 0
  else
    call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:2),dimid(i,4)/), varid(i:i) )
  end if

  ! 4D variables
  do i = 11,12
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,:), varid(i:i) )
    end if
  end do

  ! rest of 3D lon+lat+time variables
  do i = 13,nvar
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:2),dimid(i,4)/), varid(i:i) )
    end if
  end do


  ! Variables attributes:
  !----------------------

  do i = 6,nvar
    if (fnum(i)/=0) then ! if variable has to be written
      call put_att_text( fid(i), varid(i:i), (/'name'/) , varname(i:i)             )
      call put_att_text( fid(i), varid(i:i), (/'units'/), varunits(i:i)            )
      call put_att_real( fid(i), varid(i:i), missvalname(i:i) , real(missval(i:i)) ) 
      if (.not. (var_longname(i) == '')) then
          call put_att_text(fid(i), varid(i:i), (/'long_name'/), var_longname(i:i) )
      end if
    end if
  end do


  ! end of denifition :
  !--------------------

  k = 0
  do i = 5,nvar
    if (fnum(i) > k) then ! if "new" file
      call enddef(fid(i))
      k = fnum(i)
    end if
  end do


  ! put litho_frac and slope variables (#6 and #9, only time-independant variables):
  !---------------------------------------------------------------------------------

  ! Put fillvalue
  do j = 1,nlon*nlat
    if (area(j)==0) then
      litho_frac(:,j) = missval(6)
      slope(j)        = missval(9)
    end if
  end do

  i = 6
  if (fnum(i)>0) then ! if file exist
    call put_var_real3D( fid(i), varid(i), real(reshape(litho_frac, shape=(/nlon,nlat,nlit/), order=(/3,1,2/))) )
  end if

  i = 9
  if (fnum(i)>0) then ! if file exist
    call put_var_real2D( fid(i), varid(i), real(reshape(slope, shape=(/nlon,nlat/))) )
  end if


  ! close output file(s):
  !----------------------

  k = 0
  do i = 6,nvar
    if (fnum(i) > k) then ! if "new" file
      call close_file(fid(i))
      k = fnum(i)
    end if
  end do




end subroutine

end module
