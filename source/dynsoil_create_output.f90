module dynsoil_create_output_mod
implicit none

contains

subroutine dynsoil_create_output( ID,output_path,run_name, &
                                  xlevs,lon,lat,area,litho_frac,slope, fnum,fname,varname,varunits,missvalname,missval )
  use io_module, only: read_io_condition
  use netcdf_io_module, only: create, def_dim, def_var, put_att_text, put_att_real, enddef, redef, put_var_int1D, put_var_real1D, &
                              put_var_real2D, put_var_real3D, close_file
  use netcdf

  include 'coupler.inc'
  !
  integer, parameter:: nvar=25
  !
  integer, intent(in):: ID
  character(len=*), intent(in):: output_path, run_name
  double precision, intent(in), dimension(:):: xlevs, lon, lat, area, slope
  double precision, intent(in), dimension(:,:):: litho_frac
  integer, intent(out), dimension(nvar):: fnum
  character(len=*), intent(out), dimension(nvar):: fname, varname, varunits, missvalname
  double precision, intent(out), dimension(nvar):: missval
  !
  integer, dimension(nvar):: fid, varid
  !
  character(len=500), dimension(nvar):: var_longname
  integer, parameter:: maxlgth=200
  character(len=maxlgth):: dummy
  integer:: nlon, nlat, nlit, nlev, dimid(nvar,5)
  integer:: i, j, k

  ! Default name and value of missing-values (if unspecified in text input files)
  character(len=*), parameter:: defmissvalname = '_FillValue'
  double precision, parameter:: defmissval = 9.96921e+36


  !===================================== OUTPUT VARIABLES LIST: ======================================!
  !  OUTPUT: X, Y, litho, xlevs, t, area, litho_frac, slope, temp, runoff, h_soil, x_surf, tau_surf,  !
  !  var #:  1  2  3      4      5  6     7           8      9     10      11      12      13         !
  !          z, tau, Reg_prod, Reg_eros, reg_P_diss, reg_P_eros, x_surf_eros, x_P_mean, reg_mean_age, !
  !          14 15   16        17        18          19          20           21        22            !
  !          Li_Friv,  Li_Fsp,   Li_driv                                                              !
  !          23        24        25                                                                   !
  !===================================================================================================!




  ! default missing-value
  missvalname = defmissvalname
  missval = defmissval

  ! get dimension size:
  nlon = size(lon)
  nlat = size(lat)
  nlit = size(litho_frac,1)
  nlev = size(xlevs)



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%                                   getting output file names:                                   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


  ! read variables attributes (#1 to #25)
  do i=1,5 ! dimension variables
    call read_io_condition(ID, fname(i) , varname(i) , varunits(i) )
  end do
  do i=6,nvar ! rest of variables
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

  fnum(7) = 7 ! the first 6 variables are dimensions variables + area, and don't need a file
  do k = 7,nvar  ! loop on all file names
    i = 1
    do while ( fname(k) /= fname(i) ) ! loop on all previous file names
      i = i + 1
    end do
    fnum(k) = i
  end do
  do k = 1,nvar
    if (fname(k)=='') then
      ! empty file name => don't write this variable
      ! => put 0, so that it will not be created
      fnum(k) = 0
    end if
  end do


  ! Output file creation:
  !---------------------

  k = 0
  do i = 7,nvar
    if ( fnum(i) > k ) then ! if "new" file

      call create(fname(i), fid(i))
      k = fnum(i)

      ! GOBAL ATTRIBUTES
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'title'/),             (/'DynSoil module ouputs'/))
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'run_name'/),          (/trim(run_name)/))
      if (use_dynsoil_steady_state) then
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'DynSoil_mode'/),    (/'steady-state'/))
      else
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'DynSoil_mode'/),    (/'dynamic'/))
      end if
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'CO2_interpolation'/), (/interpolation_mode/))

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

  dimid = 0
  varid = 0


  ! horizontal dimensions variables:
  !---------------------------------

  k = 0
  do i = 7,nvar ! loop on all non-dimension variables nor area (as there are all horizontally-defined)
    if ( fnum(i) > k ) then
      call def_dim( fid(i), varname(1:2), (/nlon,nlat/),  dimid(i,1:2)                )
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
    elseif(fnum(i)>0) then ! if old file
      dimid(i,1:2) = dimid(fnum(i),1:2)
    end if
  end do


  ! lithology dimension variable:
  !------------------------------

  k = 0
  i = 7 ! loop on lithology-defined variables
  if ( fnum(i) > k ) then ! if "new" file
    ! dimension and variables definition:
    call def_dim( fid(i), varname(3:3), (/nlit/),     dimid(i,3:3)  )
    call def_var( fid(i), varname(3:3), (/NF90_INT/), dimid(i,3:3), varid(3:3)    )
    ! Variable attributes:
    call put_att_text( fid(i), varid(3:3), (/'name'/), varname(3:3) )
    call enddef( fid(i) )
    ! Put variable:
    call put_var_int1D( fid(i), varid(3), (/(j,j=1,nlit)/) )
    call redef( fid(i) )
    !
    do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,3) = dimid(i,3) ; end do
    !
    k = fnum(i)
  elseif(fnum(i)>0) then ! if old file
    dimid(i,3) = dimid(fnum(i),3)
  end if
  do i = 11,25 ! loop on lithology-defined variables
    if ( fnum(i) > k ) then ! if "new" file
      ! dimension and variables definition:
      call def_dim( fid(i), varname(3:3), (/nlit/),     dimid(i,3:3)  )
      call def_var( fid(i), varname(3:3), (/NF90_INT/), dimid(i,3:3), varid(3:3)    )
      ! Variable attributes:
      call put_att_text( fid(i), varid(3:3), (/'name'/), varname(3:3) )
      call enddef( fid(i) )
      ! Put variable:
      call put_var_int1D( fid(i), varid(3), (/(j,j=1,nlit)/) )
      call redef( fid(i) )
      !
      do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,3) = dimid(i,3) ; end do
      !
      k = fnum(i)
    elseif(fnum(i)>0) then ! if old file
      dimid(i,3) = dimid(fnum(i),3)
    end if
  end do


  ! vertical dimension variable:
  !-----------------------------

  k = 0
  do i = 14,15 ! loop on 5D variables
    if ( fnum(i) > k ) then ! if "new" file
      ! dimension and variables definition:
      call def_dim( fid(i), varname(4:4), (/nlev/) , dimid(i,4:4) )
      call def_var( fid(i), varname(4:4), (/NF90_FLOAT/), dimid(i,4:4), varid(4:4)    )
      ! dimension variables attributes:
      call put_att_text( fid(i), varid(4:4), (/'axis'/), (/'Z'/)        )
      call put_att_text( fid(i), varid(4:4), (/'positive'/), (/'down'/) )
      ! Variable attributes:
      call put_att_text( fid(i), varid(4:4), (/'name'/) , varname(4:4)  )
      call put_att_text( fid(i), varid(4:4), (/'units'/), varunits(4:4) )
      call enddef( fid(i) )
      ! Put variable:
      call put_var_real1D( fid(i), varid(4), real(xlevs) )
      call redef( fid(i) )
      !
      do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,4) = dimid(i,4) ; end do
      !
      k = fnum(i)
    elseif(fnum(i)>0) then ! if old file
      dimid(i,4) = dimid(fnum(i),4)
    end if
  end do


  ! time dimension and area variables:
  !-----------------------------------

  k = 0
  do i = 9,nvar ! loop on all non-dimension variables nor area and litho_frac (as there are all time dependent)
    if ( fnum(i) > k ) then
      call def_dim( fid(i), varname(5:5), (/NF90_UNLIMITED/), dimid(i,5:5)             ) ! time dim
      call def_var( fid(i), varname(5:5), (/NF90_FLOAT/),     dimid(i,5:5), varid(5:5) ) ! time
      call def_var( fid(i), varname(6:6), (/NF90_FLOAT/),     dimid(i,1:2), varid(6:6) ) ! area
      ! dimension variables attributes:
      call put_att_text( fid(i), varid(5:5), (/'axis'/), (/'T'/)                      )
      ! Variable attributes:
      call put_att_text( fid(i), varid(5:6), (/'name'/) , varname(5:6)                ) ! time and area
      call put_att_text( fid(i), varid(5:6), (/'units'/), varunits(5:6)               ) ! time and area
      call put_att_real( fid(i), varid(6:6), missvalname(6:6) , real(missval(6:6))    ) ! area
      if (.not. (var_longname(6) == '')) then
          call put_att_text(fid(i), varid(6:6), (/'long_name'/), var_longname(6:6) ) ! area
      end if
      call enddef( fid(i) )
      ! Put area variable:
      call put_var_real2D( fid(i), varid(6), real(1e12*reshape(area, shape=(/nlon,nlat/))) )
      !            area conversion: 1e6km2 -> m2  ^^^^
      call redef( fid(i) )
      !
      do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,5) = dimid(i,5) ; end do
      !
      k = fnum(i)
    elseif(fnum(i)>0) then ! if old file
        dimid(i,5:5) = dimid(fnum(i),5:5)
    end if
  end do



  !++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Other variables [need to be defined one time only] !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++!


  ! Variables definition:
  !----------------------

  ! "3D horizontal + litho" variables (litho_frac)
  if (fnum(7)==0) then ! => not writing this variable
    fid(7) = 0
    varid(7) = 0
  else
    call def_var( fid(7), varname(7:7), (/NF90_FLOAT/), dimid(7,1:3), varid(7:7) )
  end if

  ! "2D horizontal" variables (slope)
  if (fnum(8)==0) then ! => not writing this variable
    fid(8) = 0
    varid(8) = 0
  else
    call def_var( fid(8), varname(8:8), (/NF90_FLOAT/), dimid(8,1:2), varid(8:8) )
  end if

  ! "3D horizontal + time" variables (climate)
  do i = 9,10
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:2),dimid(i,5)/), varid(i:i) )
    end if
  end do

  ! "4D horizontal + litho + time" variables
  do i = 11,13
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:3),dimid(i,5)/), varid(i:i) )
    end if
  end do

  ! 5D variables ("horizontal + litho + vertical + time")
  do i = 14,15
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,:),   varid(i:i) )
    end if
  end do

  ! rest of 4D "horizontal + litho + time" variables
  do i = 16,nvar
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), (/dimid(i,1:3),dimid(i,5)/), varid(i:i) )
    end if
  end do


  ! Variables attributes:
  !----------------------

  do i = 7,nvar
    if (fnum(i)/=0) then ! if variable has to be written
      call put_att_text( fid(i), varid(i:i), (/'name'/) , varname(i:i)             )
      call put_att_text( fid(i), varid(i:i), (/'units'/), varunits(i:i)            )
      call put_att_real( fid(i), varid(i:i), missvalname(i:i) ,real( missval(i:i)) ) 
      if (.not. (var_longname(i) == '')) then
          call put_att_text(fid(i), varid(i:i), (/'long_name'/), var_longname(i:i) )
      end if
    end if
  end do


  ! end of denifition:
  !-------------------

  k = 0
  do i = 7,nvar
    if (fnum(i) > k) then ! if "new" file
      call enddef(fid(i))
      k = fnum(i)
    end if
  end do


  ! put slope and litho_frac variables (#7 and #8, only time-independant variables):
  !---------------------------------------------------------

  if (fnum(7)>0) then ! if file exist
    call put_var_real3D( fid(7), varid(7), real(reshape(litho_frac, shape=(/nlon,nlat,nlit/), order=(/3,1,2/))) )
  end if
  if (fnum(8)>0) then ! if file exist
    call put_var_real2D( fid(8), varid(8), real(reshape(slope, shape=(/nlon,nlat/))) )
  end if


  ! close output file(s):
  !----------------------

  k = 0
  do i = 7,nvar
    if (fnum(i) > k) then ! if "new" file
      call close_file(fid(i))
      k = fnum(i)
    end if
  end do




end subroutine


end module
