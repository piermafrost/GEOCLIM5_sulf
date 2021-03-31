module geoclim_write_output_mod
implicit none

contains

subroutine geoclim_write_output(t)
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, put_var_real1D
  use netcdf
  use constante, only: gotoshelf, rho_sed, PI_n_O2_atm, PI_n_CO2_atm, PI_n_rest_of_atm
  include 'combine_foam.inc'
  integer, dimension(nGEOoutvar):: nt, fileid, varid
  integer:: ierr,timevarid(1),dimid(1),fnum(nGEOoutvar)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!           offline computations:           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xPOPexport=0.
    do j0=1,nsurface
        j = jbox_surface(j0)
        xPOPexport=xPOPexport + fsink(j)*var(7,j)*vol(j)
        !totalCO2_exch=totalCO2_exch+fco2atm_ocean(j)
    enddo

    vol_tot=0.
    var_tot(:)=0.
    ph_tot=0.
    omega_tot=0.
    temp_tot=0.
    salin_tot=0.
    freef_tot=0.
    fodc_tot=0.
    fdep_tot=0.
    fco2atm_ocean_tot=0.
    F_seafloor_cdiss_tot=0.
    do j=1,nbasin-1
        vol_tot = vol_tot+vol(j)
        var_tot(1:11) = var_tot(1:11) + var(1:11,j)*vol(j)  !! concentration
        var_tot(20) = var_tot(20) + var(20,j)*vol(j)        !! of elements
        var_tot(13) = var_tot(13) + var(13,j)*var(1,j)*vol(j) ! DIC d13C
        var_tot(14) = var_tot(14) + var(14,j)*var(10,j)*vol(j) ! PIC d13C
        var_tot(15) = var_tot(15) + var(15,j)*var(9,j)*vol(j) ! POC d13C
        var_tot(17:18) = var_tot(17:18) + var(17:18,j)*var(5:6,j)*vol(j) ! diss & PIC Sr iso ration
        var_tot(19) = var_tot(19) + var(19,j)*var(20,j)*vol(j) ! d7Li
        ph_tot     = ph_tot + pH(j)*vol(j)
        omega_tot  = omega_tot + omega(j)*vol(j)
        temp_tot   = temp_tot + temp_box(j)*vol(j)
        salin_tot  = salin_tot + salin(j)*vol(j)
        freef_tot = freef_tot+freef(j)
        fodc_tot  = fodc_tot+fodc(j)
        fco2atm_ocean_tot    = fco2atm_ocean_tot + fco2atm_ocean(j)
        F_seafloor_cdiss_tot = F_seafloor_cdiss_tot + F_seafloor_cdiss(j)
    end do
    ! Carbonate deposition: sum "fsink*[PIC]" on sedimentary reservoirs
    ! + remove epicont surf fraction of fsink that did not go into shelf (1-gotoshelf)
    do j0=1,nsedi
        j = jbox_sedi(j0)
        fdep_tot = fdep_tot + fsink_inorg(j)*var(10,j)*(vol(j))
    end do
    do j0=1,nappcont-1!<- skip atmospheric box, counted in "appcont"
        j = jbox_appcont(j0)
        fdep_tot = fdep_tot - (1-gotoshelf)*fsink_inorg(j)*var(10,j)*(vol(j))
    end do

    var_tot(13) = var_tot(13) / var_tot(1) ! DIC d13C
    var_tot(14) = var_tot(14) / var_tot(10) ! PIC d13C
    var_tot(15) = var_tot(15) / var_tot(9) ! POC d13C
    var_tot(17:18) = var_tot(17:18) / var_tot(5:6) ! diss & PIC Sr iso ration
    var_tot(19) = var_tot(19) / var_tot(20) ! d7Li
    var_tot(1:11) = var_tot(1:11)/vol_tot  !! concentration
    var_tot(20) = var_tot(20)/vol_tot      !! of elements
    ph_tot    = ph_tot/vol_tot
    omega_tot = omega_tot/vol_tot
    temp_tot  = temp_tot/vol_tot
    salin_tot = salin_tot/vol_tot

    ! Present-day normalized mol of O2 and CO2
    O2_atm_level  = var(11,nbasin) / PI_n_O2_atm
    CO2_atm_level = var(12,nbasin) / PI_n_CO2_atm

    ! convert mol to volume fraction (= molar fraction) for 02 and CO2
    O2_atm_conc   =  ( var(11,nbasin) / ( PI_n_rest_of_atm + var(11,nbasin) + var(12,nbasin) ) )*100 ! (%)
    CO2_atm_conc  =  ( var(12,nbasin) / ( PI_n_rest_of_atm + var(11,nbasin) + var(12,nbasin) ) )*1e6 ! (ppm)



!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!

    fnum = GEO_ofile_num

    ! open file(s) and write time variable:
    !
    k = 0
    !
    do i = 6,nGEOoutvar ! loop on all variables
      if (fnum(i) > 0) then
        if ( fnum(i) > k ) then
          call open_file( GEO_ofile_name(i) , j , mode=NF90_WRITE )
          fileid(i) = j
          call inquire_dim( fileid(i) , GEO_varout_name(2:2) , dimid )
          ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
          nt(i) = nt(i) + 1
          call inquire_var( fileid(i) , GEO_varout_name(2:2) , timevarid )
          call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
          k = fnum(i)
        else
          fileid(i) = fileid(fnum(i))
          nt(i) = nt(fnum(i))
        end if
        call inquire_var( fileid(i) , GEO_varout_name(i:i) , varid(i:i) )
      end if
    end do



    !!!!!!!!!!!!!!!!!!!!
    ! write variables: !
    !!!!!!!!!!!!!!!!!!!!

    do i = 6,25
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(var(i-5,:))                  , (/1,nt(i)/), (/nbasin,1/)  )
    end do
    !
    i = 26
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(h2co3(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 27
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(hco3(:))                     , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 28
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(co3(:))                      , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 29
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(dh2co3(:))                   , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 30
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(dhco3(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 31
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(dco3(:))                     , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 32
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(ph(:))                       , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 33
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(omega(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 34
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(temp_box(1:nbasin-1)-273.15), (/1,nt(i)/), (/nbasin-1,1/) )
    !
    i = 35
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(salin(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 36
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(dplysc(:))                   , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 37
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(dplysa(:))                   , (/1,nt(i)/), (/nbasin,1/)  )
    !
    do i = 38,48
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(var_tot(i-37))/)           , (/nt(i)/), (/1/)           )
    end do
    !
    do i = 49,51
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(var_tot(i-36))/)           , (/nt(i)/), (/1/)           )
    end do
    !
    do i = 52,55
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(var_tot(i-35))/)           , (/nt(i)/), (/1/)           )
    end do
    !
    i = 56
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(ph_tot)/)                  , (/nt(i)/), (/1/)           )
    !
    i = 57
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(omega_tot)/)               , (/nt(i)/), (/1/)           )
    !
    i = 58
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(temp_tot-273.15)/)         , (/nt(i)/), (/1/)           )
    !
    i = 59
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(salin_tot)/)               , (/nt(i)/), (/1/)           )
    !
    i = 60
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(O2_atm_level)/)            , (/nt(i)/), (/1/)           )
    !
    i = 61
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(CO2_atm_level)/)           , (/nt(i)/), (/1/)           )
    !
    i = 62
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(O2_atm_conc)/)             , (/nt(i)/), (/1/)           )
    !
    i = 63
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(CO2_atm_conc)/)            , (/nt(i)/), (/1/)           )
    !
    i = 64
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(xPOPexport/1.d+12)/)       , (/nt(i)/), (/1/)           )
    !
    i = 65
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fanthros/1.d+12)/)         , (/nt(i)/), (/1/)           )
    !
    i = 66
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fco2atm_ocean(:))            , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 67
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fco2atm_ocean_tot)/)       , (/nt(i)/), (/1/)           )
    !
    i = 68
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(finorgC(:))                  , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 69
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fdissol_carb(3)*var(10,3)*vol(3))/) &
                                                                                                      , (/nt(i)/), (/1/)           )
    !
    i = 70
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fsilw+fbasw)/)             , (/nt(i)/), (/1/)           )
    !
    i = 71
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fbasw)/)                   , (/nt(i)/), (/1/)           )
    !
    i = 72
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fcarbw)/)                  , (/nt(i)/), (/1/)           )
    !
    i = 73
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fkerw)/)                   , (/nt(i)/), (/1/)           )
    !
    i = 74
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(freef(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 75
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(freef_tot)/)               , (/nt(i)/), (/1/)           )
    !
    i = 76
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fdep_tot)/)                , (/nt(i)/), (/1/)           )
    !
    i = 77
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fodc(:))                     , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 78
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fodc_tot)/)                , (/nt(i)/), (/1/)           )
    !
    i = 79
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fsink_inorg(:)*var(10,:)*vol(1:nbasin)) &
                                                                                                      , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 80
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fpw)/)                     , (/nt(i)/), (/1/)           )
    !
    i = 81
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fbioC(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 82
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(F_seafloor_cdiss(:))         , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 83
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(F_seafloor_cdiss_tot)/)    , (/nt(i)/), (/1/)           )
    !
    i = 84
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(ftrap)/)                   , (/nt(i)/), (/1/)           )
    !
    i = 85
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fCO2_crust(:))               , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 86
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fSO4_basin(:))               , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 87
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fSO4_crust(:))               , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 88
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(FrivLi)/)                  , (/nt(i)/), (/1/)           )
    !
    i = 89
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(dLiriv)/)                  , (/nt(i)/), (/1/)           )
    !
    i = 90
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(total_cont_POC_export)/)   , (/nt(i)/), (/1/)           )
    !
    i = 91
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fvol+sum(fmor)+ftrap+fanthros)/)  , (/nt(i)/), (/1/)    )
    !
    i = 92
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(discharge)/)               , (/nt(i)/), (/1/)           )
    !
    i = 93
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(tss)/)                     , (/nt(i)/), (/1/)           )
    !
    i = 94
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(epsiC(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 95
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(epsiCont)/)                , (/nt(i)/), (/1/)           )
    !
    i = 96
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(temp_box(nbasin:nbasin))     , (/nt(i)/), (/1/)           )
    !
    i = 97
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(ws(:))                       , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 98
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(rho_sed*fin_sed(:))          , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 99
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(Corg_BE(:))                  , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 100
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fodp(:))                     , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 101
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fphos(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 102
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fhydP(:))                    , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 103
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fcarbsulfw+fsilsulfw+fH2SO4sulfw)/) , (/nt(i)/), (/1/)  )
    !
    i = 104
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fsilsulfw)/)               , (/nt(i)/), (/1/)           )
    !
    i = 105
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  (/real(fcarbsulfw)/)              , (/nt(i)/), (/1/)           )
    !
    i = 106
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fSulfRed(:))                 , (/1,nt(i)/), (/nbasin,1/)  )
    !
    i = 107
      if (fnum(i)>0)    call put_var_real1D(  fileid(i), varid(i),  real(fO2_odc(:))                  , (/1,nt(i)/), (/nbasin,1/)  )

  !======================!
  ! Output file closing: !
  !======================!

  k = 0
  do i = 6,nGEOoutvar
    if ( fnum(i) > k ) then
      call close_file(fileid(i))
      k = fnum(i)
    end if
  end do



end subroutine


end module
