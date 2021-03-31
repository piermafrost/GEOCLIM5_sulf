module geoclim_create_output_mod
implicit none

contains

subroutine geoclim_create_output(ID, output_path, run_name, vol, surf, surf_sedi, fname, fnum, varname, missval)
  use io_module, only: read_io_condition
  use netcdf_io_module, only: create, def_dim, def_var, put_att_text, put_att_real, enddef, redef, put_var_real1D, &
                              put_var_int1D, close_file
  use netcdf
  include 'coupler.inc'
  integer, parameter:: nvar=107
  integer, intent(in):: ID
  character(len=*), intent(in):: output_path, run_name
  double precision, dimension(:), intent(in):: vol, surf, surf_sedi
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
  integer:: nbasin, nlon, nlat, dimid(nvar,2)
  integer:: i, j, k

  ! Default name and value of missing-values (if unspecified in text input files)
  character(len=*), parameter:: defmissvalname = '_FillValue'
  double precision, parameter:: defmissval = 9.96921e+36

  !===================================== OUTPUT VARIABLES LIST: =====================================!
  !  OUTPUT: box, time, vol, surf, surf_sed, DIC, alk, PO4, Ca, Sr, SrPIC, POP, PIP, POC, PIC, O2,   !
  !  var #:  1    2     3    4     5         6    7    8    9   10  11     12   13   14   15   16    !
  !          PCO2,    DICd13C, PICd13C, POCd13C, CO2d13C, SrISO, PICSrISO, d7Li, Li, H2CO3, HCO3,    !
  !          17       18       19       20       21       22     23        24    25  26     27       !
  !          CO3, H2CO3d13C, HCO3d13C, CO3d13C, pH, omega, temp, salin, dplysc, dplysa, DIC_G,       !
  !          28   29         30        31       32  33     34    35     36      37      38           !
  !          alk_G, PO4_G, Ca_G, Sr_G, SrPIC_G, POP_G, PIP_G, POC_G, PIC_G, O2_G, DICd13C_G,         !
  !          39     40     41    42    43       44     45     46     47     48    49                 !
  !          PICd13C_G, POCd13C_G, SrIOS_G, PICSrISO_G, d7Li_G, Li_G, pH_G, omega_G, temp_G,         !
  !          50         51         52       53          54      55    56    57       58              !
  !          salin_G, O2_LEV, CO2_LEV, O2_CON, CO2_CON, xPOPex, fCO2anth, FCO2AO, FCO2AO_tot,        !
  !          59       60      61       62      63       64      65        66      67                 !
  !          finorgC, Fdisscarb, Fsilw, Fbasw, Fcarbw, Fkerw, freef, freef_tot, fdep_tot, Fodc,      !
  !          68       69         70     71     72      73     74     75         76        77         !
  !          Fodc_tot, Fsinkinorg, FPw, fbioC, FseaflrCdiss, FseaflrCdiss_tot, ftrap, fCO2crust,     !
  !          78        79          80   81     82            83                84     85             !
  !          FSO4basin, FSO4crust, LiFriv, Lidriv, total_cont_POC_export, degassing, discharge, TSS, !
  !          86         87         88      89      90                     91         92         93   !
  !          epsi_C, epsi_C_cont, GMST, sed_rate, sed_flux, Corg_BE, fodp, fphos, fhydP, fPyrWth,    !
  !          94      95           96    97        98        99       100   101    102    103         !
  !          fSilSulfW, fCarbSulfW, fSulfRed, fO2_odc                                                !
  !          104        105         106       107                                                    !
  !==================================================================================================!

  integer, dimension(nvar), parameter:: indice_boxdef = &
          (/ 0,0,1,1,1, 1,1,1,1,1,   1,1,1,1,1, 1,1,1,1,1,   1,1,1,1,1, 1,1,1,1,1,   1,1,1,1,1, 1,1,0,0,0, & !#40
             0,0,0,0,0, 0,0,0,0,0,   0,0,0,0,0, 0,0,0,0,0,   0,0,0,0,0, 1,0,1,0,0,   0,0,0,1,0, 0,1,0,1,0, & !#80
             1,1,0,0,1, 1,1,0,0,0,   0,0,0,1,0, 0,1,1,1,1,   1,1,0,0,0, 1,1                                /)




  ! default missing-value
  missvalname = defmissvalname
  missval = defmissval

  ! get number of basin
  nbasin = size(surf,1)


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%                                   getting output file names:                                   %!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


  ! read variables attributes (#1 to #nvar)
  do i=1,5
    call read_io_condition(ID, fname(i) , varname(i) , varunits(i) )
  end do
  do i=6,nvar
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

  fnum(6) = 6 ! the first 5 variables are dimensions variables + volume + surface, and don't need a file
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
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'title'/),                       (/'COMBINE module ouputs'/))
      call put_att_text(fid(i), (/NF90_GLOBAL/), (/'run_name'/),                    (/trim(run_name)/))
      if (coupling_ecogeo) then
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'marine_ecology_module'/),     (/'on'/))
      else
        call put_att_text(fid(i), (/NF90_GLOBAL/), (/'marine_ecology_module'/),     (/'off'/))
      endif
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


  ! box dimension and variable
  !---------------------------

  k = 0
  do i = 6,nvar ! loop on box-defined variables (except axes, volume and surfaces)
    if (indice_boxdef(i)==1) then

      if ( fnum(i) > k ) then ! if "new" file

        call def_dim( fid(i), varname(1:1), (/nbasin/) , dimid(i,1:1) )
        call def_var( fid(i), varname(1:1), (/NF90_INT/), dimid(i,1:1), varid(1:1)    )
        ! dimension variables attributes:
        call put_att_text( fid(i), varid(1:1), (/'axis'/), (/'Y'/)        )
        ! Variable attributes:
        call put_att_text( fid(i), varid(1:1), (/'name'/) , varname(1:1)  )
        call put_att_text( fid(i), varid(1:1), (/'units'/), varunits(1:1) )
        call put_att_text( fid(i), varid(1:1), (/'box_1'/),  (/'polar N surf'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_2'/),  (/'polar N deep'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_3'/),  (/'mid-lat surf'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_4'/),  (/'thermocline'/)  )
        call put_att_text( fid(i), varid(1:1), (/'box_5'/),  (/'mid-lat deep'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_6'/),  (/'epicont surf'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_7'/),  (/'epicont deep'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_8'/),  (/'polar S surf'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_9'/),  (/'polar S deep'/) )
        call put_att_text( fid(i), varid(1:1), (/'box_10'/), (/'atmosphere'/)   )

        do j = fnum(i),i-1 ; if (fnum(j)==fnum(i)) dimid(j,1) = dimid(i,1) ; end do

        k = fnum(i)

      elseif(fnum(i)/=0) then ! if "old" file
        dimid(i,1) = dimid(fnum(i),1)
      end if


    end if
  end do


  ! time dimension and variable:
  !-----------------------------

  k = 0
  do i = 6,nvar ! loop on all variables (except axes and volume)

    if ( fnum(i) > k ) then ! if "new" file

      ! time dimension (all variables are time-defined):
      call def_dim( fid(i), varname(2:2), (/NF90_UNLIMITED/) , dimid(i,2:2) )
      call def_var( fid(i), varname(2:2), (/NF90_FLOAT/), dimid(i,2:2), varid(2:2)    )
      ! dimension variables attributes:
      call put_att_text( fid(i), varid(2:2), (/'axis'/), (/'T'/)        )
      ! Variable attributes:
      call put_att_text( fid(i), varid(2:2), (/'name'/) , varname(2:2)  )
      call put_att_text( fid(i), varid(2:2), (/'units'/), varunits(2:2) )
    
      k = fnum(i)

    elseif(fnum(i)/=0) then ! if "old" file
      dimid(i,2) = dimid(fnum(i),2)
    end if

  end do


  ! volume and surface variable:
  !-----------------------------

  k = 0
  do i = 6,nvar ! loop on box-defined variables (except axes and volume)
    if (indice_boxdef(i)==1) then

      if ( fnum(i) > k ) then ! if "new" file

        call def_var( fid(i), varname(3:3), (/NF90_FLOAT/), dimid(i,1:1), varid(3:3)    )
        call put_att_text( fid(i), varid(3:3), (/'name'/) , varname(3:3)  )
        call put_att_text( fid(i), varid(3:3), (/'units'/), varunits(3:3) )

        call def_var( fid(i), varname(4:4), (/NF90_FLOAT/), dimid(i,1:1), varid(4:4)    )
        call put_att_text( fid(i), varid(4:4), (/'name'/) , varname(4:4)  )
        call put_att_text( fid(i), varid(4:4), (/'units'/), varunits(4:4) )

        call def_var( fid(i), varname(5:5), (/NF90_FLOAT/), dimid(i,1:1), varid(5:5)    )
        call put_att_text( fid(i), varid(5:5), (/'name'/) , varname(5:5)  )
        call put_att_text( fid(i), varid(5:5), (/'units'/), varunits(5:5) )

        k = fnum(i)

      end if

    end if
  end do




  !++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Other variables [need to be defined one time only] !
  !++++++++++++++++++++++++++++++++++++++++++++++++++++!


  ! Variables definition:
  !----------------------

  do i = 6,nvar
    if (fnum(i)==0) then ! => not writing this variable
      fid(i) = 0
      varid(i) = 0
    else
      if (indice_boxdef(i)==0) then ! if box-defined variable:
        call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,2:2), varid(i:i) )
      else
        call def_var( fid(i), varname(i:i), (/NF90_FLOAT/), dimid(i,1:2), varid(i:i) )
      end if
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
  do i = 6,nvar
    if (fnum(i) > k) then ! if "new" file
      call enddef(fid(i))
      k = fnum(i)
    end if
  end do


  ! put box, volume and surface variables:
  !---------------------------------------

  k = 0
  do i = 6,nvar
    if (indice_boxdef(i)==1) then
      if ( fnum(i) > k ) then ! if "new" file
        call put_var_int1D( fid(i), varid(1), (/(j,j=1,nbasin)/) )
        call put_var_real1D( fid(i), varid(3), real(vol(1:nbasin)) )
        call put_var_real1D( fid(i), varid(4), real(surf) )
        call put_var_real1D( fid(i), varid(5), real(surf_sedi) )
        k = fnum(i)
      end if
    end if
  end do


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
