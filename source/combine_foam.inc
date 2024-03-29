!   Technical parameters and variables:
!   ==================================================================

    include 'path.inc' ! => get the variable 'geoclim_path'
    character(len=100)    :: run_name, check_run_name, killing_file_name
    character(len=500)    :: output_path
    character(len=30)     :: land_input_mode
    integer, dimension(5) :: ERROR_HANDLING_OPTION

    include 'coupler.inc' ! coupling parameters: 'coupling_dynsoil' and 'coupling_ecogeo' plus 'interpolation_mode' 

    include 'shape.inc'
    ! => define: number of boxes ('nbasin')
    !            number of CO2 levels ('nclimber')
    !            geographic dimensions ('nlon' and 'nlat') -> number of pixels
    !            number of lithology classes ('nlitho'),
    !            number of DynSoil levels ('nDSlev') 
    integer, parameter :: nvar_real=20
    integer, parameter :: npixel=nlon*nlat
    integer, parameter :: nveg=80    !number of primary producers, including empty slots 80 60
    integer, parameter :: npred=10   !number of predator type II  10 15
    integer, parameter :: nspred=5    !number of p1edator type III 5 10
    integer, parameter :: nsspred=3  !number of predator type IV  3
    integer, parameter :: nclado_veg=0
    integer, parameter :: nequat=nveg+npred+nspred+nsspred  !number of equation for the biodiversity model
    integer, parameter :: nclass=4   !obsolete total number of bio types
    integer, parameter :: nprint=100000

    integer, parameter :: nvar=nvar_real*nbasin,ndegas=1000000,nsres=24
    integer, parameter :: nGEOoutvar=107, nBIOoutvar=7, nGEOGoutvar=15

    integer, parameter :: nDSvar = 25       ! number of dynsoil output variables
    integer, parameter :: n_peaks = 5





!   Variables:
!   ==================================================================
    logical decide,check
    integer i,j,j0,k,kl,icount,n,ishelfal,k1,k2,ncount,icompteur,kbasin
    integer i_rnd,iextinct,jmute,nrandom
    integer clock
    integer extrinsic_flag,ipeak
    integer:: convert2ascii
    integer ijump_cont_weath,icount_cont_weath, ijump_geogprint,icount_geogprint
    integer, dimension(:), allocatable :: seed
    integer, dimension(8):: computer_time
    
    integer kgrazed(1:npred+nspred+nsspred,1:nequat)
    integer kgrazedFlag(1:npred+nspred+nsspred,1:nveg+npred+nspred+nsspred,1:nbasin)
    integer f1bio(1:npred+nspred+nsspred)
    integer extinct(1:nequat,1:nbasin),iextext(1:nequat,1:nbasin)
    integer index_carb(1:nveg,1:nbasin)
    integer web(1:nbasin,1:nveg,1:nveg),iclado_veg(1:nbasin)
    integer num_pix_basin(1:nbasin),stop_clado(1:nbasin)
    integer fauna,iweb
    integer specount(1:nbasin)

    double precision ybio(1:nequat),ybio_0(1:nequat),dydxbio(1:nequat),youtbio(1:nequat)
    double precision ybio_steady(1:nequat,1:nbasin)
    double precision yprint(1:nequat),percent(1:nequat)
    double precision alphabio(1:nequat),omegabio(1:nequat),birth(1:nequat,1:nbasin)
    double precision sigma2(1:nequat,1:nbasin),death(1:nequat,1:nbasin),temp0(1:nequat,1:nbasin)
    double precision sigma3(1:nequat,1:nbasin),pH0(1:nequat,1:nbasin)
    double precision r_ar(10)
    double precision grazing_efficiency(1:npred+nspred+nsspred,1:nequat,1:nbasin)
    double precision grazing(1:npred+nspred+nsspred,1:nequat)
    double precision alphabio_prim(1:nbasin),ratio_survive(1:nequat)
    double precision limiting_factor(1:nbasin)
    
    double precision rbio,xbio,biotemp,biopH,furious,revolution,rshift,rextext
    double precision a,b,c,a2,b2,c2
    double precision sigma,sigma_p
    double precision biomass_veg,biomass_carb_prod
    double precision total_bio(1:nclass),grazing_tot
    double precision var_bio(1:nequat,1:nbasin)
    double precision fbiodivC(1:nbasin)
    double precision temp_adapt(1:nbasin),pH_adapt(1:nbasin)
    double precision carb_ratio(1:nbasin),interact(1:nbasin,1:nveg,1:nveg)
    double precision interact_ave(1:nbasin,1:nveg)

    double precision revol,revol_type

    double precision tstart_deg2,tend_deg2
    double precision out,O2,shelfal,xkill,sluggish
    double precision DOA(1:nbasin),cp_burial(1:nbasin),t
    double precision fCO2atm_ocean(1:nbasin),po2,fC13atm_ocean(1:nbasin), &
                     fC13ocean_atm(1:nbasin)
    double precision oce_surf(1:nbasin), oce_surf_tot, temp_box(1:nbasin), press_box(1:nbasin), surf_sedi(1:nbasin), &
                     sedim_fract(nbasin), sedim_capacity(nbasin)
    integer, dimension(nbasin) :: indice_deep, indice_sedi, app_cont,indice_thermo, indice_surface, indice_epicont, indice_polar
    integer :: ndeep, nnodeep, nsedi, nnosedi, nthermo, nnothermo, nsurface, nnosurface, &
               nepicont, nnoepicont, npolar, nnopolar, nappcont, nnoappcont, nsediepicontsurf, &
               nsediepicontnosurf, nsedinoepicont, nsurfappcont, nsurfnoappcont, nnosurfappcont, &
               nnosurfnoappcont, nnosurfbelappcont, nnosurfnobelappcont
    integer, dimension(nbasin) :: jbox_deep, jbox_nodeep, jbox_sedi, jbox_nosedi, jbox_thermo, jbox_nothermo, &
                                  jbox_surface, jbox_nosurface, jbox_epicont, jbox_noepicont, jbox_polar, jbox_nopolar, &
                                  jbox_appcont, jbox_noappcont, jbox_sediepicontsurf, jbox_sediepicontnosurf, &
                                  jbox_sedinoepicont, jbox_surfappcont, jbox_surfnoappcont, jbox_nosurfappcont, &
                                  jbox_nosurfnoappcont, jbox_nosurfbelappcont, jbox_nosurfnobelappcont
    double precision fbioP(1:nbasin),fbioC(1:nbasin),finorgC(1:nbasin), &
                     finorgP(1:nbasin)
    double precision fdownt,fdowntP
    double precision dc13PTAS,dc13PTAD,dc13PTACS,dc13PTACD,dc13PTAT &
                     ,dco3(1:nbasin),dhco3(1:nbasin),dh2co3(1:nbasin) &
                     ,dc13atm,dctrap
    double precision epsiC(1:nbasin)
    double precision epsiCont
    double precision Toceclimber(1:nbasin,1:nclimber)
    double precision Tairclimber(1:npixel,1:nclimber),Runclimber(1:npixel,1:nclimber),GMSTclimber(1:nclimber)
    double precision areaclimber(1:npixel),co2climber(1:nclimber),areatot
    double precision areaEarth(1:npixel),area_ocean(1:nbasin)
    double precision hdeepC
    double precision fvol,fmor(1:nbasin),ftrap,fSO4_atmos,fSO4_ocean,fSO4_deep,fSO4_basin(1:nbasin) &
                     ,fSO4_crust(1:nbasin),fCO2_atmos,fCO2_crust(1:nbasin),fCO2_deep
    double precision clo,phosss
    double precision var(1:nvar_real,1:nbasin)
    double precision volin,xMORin,xlumin,berner,fvegin,oblin,shells &
                     ,xnoorg,xdrift,xsolar,temcondi,async,ageYprint &
                     ,oxy_acc_fact,sulf_acc_fact,isolver
    double precision fanthros
    double precision areaPTAS,areaPTACS,dcdiss,aconv,fevapnow &
                     ,fog,ak1soil,ak2soil,akspsoil,betasoil
    double precision freef(1:nbasin),fodc(1:nbasin),fodp(1:nbasin) &
                     ,freefP(1:nbasin),f_diss_c(1:nbasin),f_diss_a(1:nbasin) &
                     ,fdepa(1:nbasin),fdepc(1:nbasin),fcarbP_dep(1:nbasin),fcarb_dep(1:nbasin) &
                     ,fdissol_carb(1:nbasin),fdissol_carbP(1:nbasin) &
                     ,ws(nbasin),fin_sed(nbasin),Corg_BE(nbasin)
    double precision F_seafloor_cdiss(1:nbasin)
    double precision tend  
    !double precision deephyp(13),shyp(13) ! obsolete
    double precision dplysc(1:nbasin),dplysa(1:nbasin)
    double precision value_keep1(1:npixel),value_keep2(1:npixel)
    double precision fphos(1:nbasin)
    double precision fhydP(1:nbasin)
    double precision pco2soil,Hsoil,Pmax
    double precision rC_Corg(1:nbasin),reff(1:nbasin)
    double precision roxyd(1:nbasin),ranoxyd(1:nbasin)
    double precision runclim(1:npixel), discharge
    double precision tss
    double precision rsw
    double precision Tclim(1:npixel)
    double precision salin(1:nbasin)
    double precision omega(1:nbasin),omega_ara(1:nbasin)
    double precision St,S
    double precision h2co3(1:nbasin),hco3(1:nbasin),co3(1:nbasin), &
                     ch(1:nbasin),ak1(1:nbasin),ak2(1:nbasin), &
                     akb(1:nbasin),akc(1:nbasin),bco2(1:nbasin),ph(1:nbasin) &
                     ,pco2_dissous(1:nbasin)
    double precision ts,tsta,tfi,xjump,ifog
    double precision rSrdep(1:nbasin)
    double precision fsink(1:nbasin),fsink_inorg(1:nbasin)
    double precision SO4_conc(1:npixel),SO4soil
    double precision fodt,ept
    double precision tdata_deg(1:ndegas),data_deg(1:ndegas)
    double precision fracveg
    double precision vol(1:nvar)
    double precision fsilw,fcarbw,fkerw,fpw,fbasw,fpyrw,ftrapw,fcarbsulfw,fsilsulfw,fH2SO4sulfw,faddsulfw ! add carbonate and silicate weathering by sulphuric acid
    ! Sulfate Reduction:
    double precision:: fodc_noSR_m2, fSulfRed_m2, fO2_odc_m2
    double precision, dimension(nbasin):: fSulfRed, fO2_odc
    double precision total_cont_POC_export,reg_eros_galy_unit
    double precision F(1:nbasin,1:nbasin)
    double precision R(1:nvar_real,1:nbasin)
    double precision ATe(1:npixel),BTe(1:npixel)
    double precision ARu(1:npixel),BRu(1:npixel)
    double precision litho_frac(1:nlitho,1:npixel)
    double precision weighted_rsw
    double precision ref_x_axis(nlon), ref_y_axis(nlat)



    double precision y(1:nvar),y_0(1:nvar),dydx(1:nvar),yout(1:nvar)
    double precision Corg_hml(1:nbasin),Corg_hsr(1:nbasin), &
                     fodcm2(1:nbasin),fin(1:nbasin),fin0,fsedim,fexport
    double precision area_cal(1:nbasin),area_ara(1:nbasin)
    double precision fsil(1:npixel),fcarb(1:npixel),fker(1:npixel)
    double precision wth_allsil(1:npixel),wth_litho_wgh(1:nlitho,1:npixel),wth_litho(1:nlitho,1:npixel)
    double precision POC_export_rate(1:npixel)
    double precision fp(1:npixel),fcarbw_pla(1:npixel)
    double precision areagra(1:npixel),fbas(1:npixel)
    double precision Runinterd(1:nclimber),Tinterd(1:nclimber)
    double precision edb(1:nbasin),ebc(1:nbasin)
    double precision beta(1:nbasin),diff(1:nvar_real,1:nbasin),bilan_water(1:nbasin)

    double precision runoff,Patm,temperature,temper
    double precision x1,x2,xacc,xn,fx,dfx
    double precision h2co3soil,hco3soil,co3soil,calcium,xmass_carbonate
    
    double precision p,totarea,pco2,fco2,fe
    double precision facgra,facbas,facker,carb_weath_conc
    double precision closed,x,xnum,xden,somme,xlogK,dum,RT0
    double precision pco2pal,fco2v,fev,fco2nv,fenv,time,pmin,prmax,aco2
    double precision htot,areal,dummy,eps,stept,compteur,tbegin,htry
    double precision sum
    double precision salsoil,prsoil,akbsoil,akcsoil,aksp,xo,fsoil,dfsoil
    double precision pr,fph,xv,dfph
    double precision xPOPexport,totalCO2_exch,freef_tot,fodc_tot,fdep_tot,fco2atm_ocean_tot,F_seafloor_cdiss_tot
    double precision var_tot(nvar_real),ph_tot,vol_tot,omega_tot,temp_tot,salin_tot,&
                     O2_atm_level,CO2_atm_level,O2_atm_conc,CO2_atm_conc
    
    double precision tka
    double precision ftrap1a,ftrap1b,ftrap2,ftrap3,ftrap4

    ! Lithium cycle
    double precision FrivLi,dLiriv,Li_dep
    double precision FhydLi(1:nbasin)

    character(len=500), dimension(nGEOoutvar) :: GEO_ofile_name, GEO_varout_name
    character(len=500), dimension(nBIOoutvar) :: BIO_ofile_name, BIO_varout_name
    character(len=500), dimension(nGEOGoutvar) :: GEOG_ofile_name, GEOG_varout_name
    integer :: GEO_ofile_num(nGEOoutvar), BIO_ofile_num(nBIOoutvar), GEOG_ofile_num(nGEOGoutvar)
    double precision :: GEO_varout_missval(nGEOoutvar), BIO_varout_missval(nBIOoutvar), GEOG_varout_missval(nGEOGoutvar)
    character(len=500):: dummychar, filename
    character(len=100):: varname, fillvalname

    ! continental variables
    integer, dimension(npixel):: list_cont_pixel
    integer:: ncontpxl

    ! dynsoil module variables:
    character(len=30):: DS_init_mode
    integer :: ijump_DS_integration, ijump_DS_print, icount_DS_int, icount_DS_pri
    double precision :: DS_timestep, xlevs(nDSlev)
    double precision, dimension(nlitho,npixel) :: reg_thick,reg_x_surf,reg_tau_surf,reg_P_vol,reg_mean_age,reg_x_mean
    double precision, dimension(nlitho,npixel) :: reg_thick_0,reg_x_surf_0,reg_tau_surf_0,reg_P_vol_0
    double precision, dimension(nlitho,npixel) :: reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros
    double precision, dimension(nlitho,npixel) :: reg_Li_Friv, reg_Li_Fsp, reg_Li_driv
    double precision, dimension(nDSlev,nlitho,npixel) :: reg_z_prof,reg_tau_prof
    double precision, dimension(nDSlev,nlitho,npixel) :: reg_z_prof_0,reg_tau_prof_0
    double precision, dimension(npixel) :: slope, reg_eros_lithmean
    logical, dimension(npixel) :: DS_missingpoints
    integer, dimension(nlitho,npixel) :: reg_ktop, reg_ktop_0
    integer, dimension(nDSvar) :: DS_ofile_num
    character(len=500), dimension(nDSvar) :: DS_ofile_name
    character(len=100), dimension(nDSvar) :: DS_varout_name
    character(len=50), dimension(nDSvar) :: DS_varout_units, DS_varout_missvalname
    character(len=200) :: DS_restart_name
    double precision, dimension(nDSvar) :: DS_varout_missval
    double precision :: DS_slopemissval
    
    
    !LIP degassing
    double precision, dimension(n_peaks) :: pulse_start,pulse_end
    double precision, dimension(n_peaks) :: dc13_peak,amount_peak
    double precision, dimension(n_peaks) :: peak_duration

    ! vegetation variables:
    integer :: coupling_veget, icount_veget, ijump_veget, veget_occup(npixel)
    double precision :: veget_factor(npixel), veget_eros_factor(npixel)

    
    data temp_adapt /288.15,275.15,293.15,999,999,293.15,999,275.15,999,999/
    data pH_adapt /7.949,7.50,8.076,999,999,8.015,999,7.955,999,999/



!   time constraint for CAMP degassing
    data pulse_start    /0, 175e3, 195e3, 265e3, 623e3/  !yrs
    data pulse_end      /50e3, 195e3, 265e3, 285e3, 643e3/  !yrs
    data dc13_peak      /-13e-3, -41e-3, -5e-3, -41e-3, -5e-3/
    data amount_peak    /8800, 4800, 4000, 2500, 4000/  !units=Gt/C
    

!   Global variables:
!   ==================================================================    

    common /IOcond/ run_name, output_path, killing_file_name, land_input_mode
    common /error/ ERROR_HANDLING_OPTION
    common /convert/ convert2ascii
    common /jump/ ijump_cont_weath,icount_cont_weath, ijump_geogprint,icount_geogprint
    common /anoxia/ DOA,cp_burial
    common /axis/ ref_x_axis, ref_y_axis
    common /atmosphere/ fCO2atm_ocean,po2,fC13atm_ocean,fC13ocean_atm
    common /basin_parameter/ oce_surf, indice_deep,indice_sedi, surf_sedi, sedim_fract, sedim_capacity, &
                             temp_box, press_box,app_cont,indice_thermo, &
                             indice_surface,indice_epicont, &
                             indice_polar,oce_surf_tot
    common /basin_ndex/ ndeep, nnodeep, nsedi, nnosedi, nthermo, nnothermo, nsurface, nnosurface, &
                        nepicont, nnoepicont, npolar, nnopolar, nappcont, nnoappcont, nsediepicontsurf, &
                        nsediepicontnosurf, nsedinoepicont, nsurfappcont, nsurfnoappcont, nnosurfappcont, &
                        nnosurfnoappcont, nnosurfbelappcont, nnosurfnobelappcont
    common /basin_jbox/ jbox_deep, jbox_nodeep, jbox_sedi, jbox_nosedi, jbox_thermo, jbox_nothermo, &
                        jbox_surface, jbox_nosurface, jbox_epicont, jbox_noepicont, jbox_polar, jbox_nopolar, &
                        jbox_appcont, jbox_noappcont, jbox_sediepicontsurf, jbox_sediepicontnosurf, jbox_sedinoepicont, &
                        jbox_surfappcont, jbox_surfnoappcont, jbox_nosurfbelappcont, jbox_nosurfnobelappcont
    common /biopump/ fbioP,fbioC,finorgC,finorgP
    common /biomass1/ ybio_steady,extinct,web,iextext
    common /biomass2/ carb_ratio,index_carb,interact,interact_ave
    common /bio_div1/ alphabio,omegabio,sigma2,temp0,birth,death
    common /bio_div2/ sigma3,pH0,fauna,iclado_veg
    common /calcite_diss_therm/ fdownt,fdowntP
    common /carb13/ dco3,dhco3,dh2co3 &
                    ,dc13atm,dctrap,edb
    common /Cfrac/ epsiCont,epsiC
    common /CLIMBER1/ Toceclimber
    common /CLIMBER2/ Tairclimber,Runclimber,GMSTclimber
    common /CLIMBER3/ areaclimber,co2climber,areaEarth,areatot
    common /coastal_depth/ hdeepC
    common /CO2_degas/ fvol,fmor,ftrap,fSO4_atmos,fSO4_ocean,fSO4_deep,fSO4_basin &
                       ,fSO4_crust,fCO2_atmos,fCO2_crust,fCO2_deep &
                       ,tstart_deg2,tend_deg2,ipeak
    common /closed_system/ clo,phosss
    common /concent/ var
    common /bio_ab/ var_bio
    common /biogeo_couple/ alphabio_prim,specount
    common /shell_alone/ shells,oxy_acc_fact,sulf_acc_fact
    common /condi/ volin,xMORin,xlumin,berner,fvegin,oblin &
                   ,xnoorg,xdrift,xsolar,temcondi,async,ageYprint &
                   ,isolver,iweb,extrinsic_flag
    common /CO2_em/ fanthros
    common /terrestrial_POC/ total_cont_POC_export, POC_export_rate
    common /cste/ areaPTAS,areaPTACS,dcdiss,aconv,fevapnow,sigma &
                  ,fog,ak1soil,ak2soil,akspsoil,betasoil
    common /deposition/ freef,fodc,fodp &
                        ,freefP,f_diss_c,f_diss_a &
                        ,fdepa,fdepc,fcarbP_dep,fcarb_dep &
                        ,fdissol_carb,fdissol_carbP,ws,fin_sed,Corg_BE
    common /dissolution/ F_seafloor_cdiss
    common /final_t/ tend  
    !common /hypsometry/ deephyp,shyp ! obsolete
    common /interclimate/ ATe,BTe,ARu,BRu
    common /litho/ litho_frac
    common /lithium/ FrivLi,dLiriv,FhydLi
    common /lysodepth/ dplysc,dplysa
    common /xnewton_keep/ value_keep1,value_keep2
    common /phosphor/ fphos
    common /Phydrot/ fhydP
    common /pressure/ pco2soil,Hsoil,Pmax
    common /ratioC_Corg/ rC_Corg,reff
    common /recycling/ roxyd,ranoxyd
    common /runoffCLIMBER/ runclim,discharge
    common /erosion/ tss
    common /r8786/ rsw
    common /temperatureCLIMBER/ Tclim
    common /saline/ salin
    common /saturation/ omega,omega_ara
    common /sedirate/ St,S
    common /species/ h2co3,hco3,co3, &
                     ch,ak1,ak2, &
                     akb,akc,bco2,ph &
                     ,pco2_dissous
    common /solver_condi/ ts,tsta,tfi,xjump,ifog
    common /str_no_cste/ rSrdep,weighted_rsw
    common /sinking/ fsink,fsink_inorg
    common /sulfate/ SO4_conc,SO4soil
    common /terorg/ fodt,ept
    common /trap/ tdata_deg,data_deg
    common /vegetal_frac/ fracveg
    common /volumes/ vol
    common /weathfl/ fsilw,fcarbw,fkerw,fpw,fbasw,fcarbsulfw,fsilsulfw,fH2SO4sulfw ! add carbonate and silicate weathering by sulphuric acid
    common /weathrt/ fsil, fbas, wth_allsil, wth_litho_wgh, wth_litho, fker, fp
    common /sulfred/ fSulfRed, fO2_odc ! variables for sulfate reduction
    common /water/ F
    common /creation/ R
    common /predator_glob1/ grazing_efficiency,grazing
    common /predator_glob2/ kgrazedFlag,i_rnd,kgrazed,iextinct
    common /shellsbug/ ishelfal
    common /outputs1/ GEO_ofile_name, GEO_varout_name
    common /outputs2/ BIO_ofile_name, BIO_varout_name
    common /outputs3/ GEOG_ofile_name, GEOG_varout_name
    common /outputs4/ GEO_ofile_num, BIO_ofile_num, GEOG_ofile_num
    common /outputs5/ GEO_varout_missval, BIO_varout_missval, GEOG_varout_missval
    common /contpxl/ list_cont_pixel, ncontpxl
    common /dynsoil_char/     DS_init_mode
    common /dynsoil_intg0D/   ijump_DS_integration, ijump_DS_print, icount_DS_int, icount_DS_pri
    common /dynsoil_dim/      DS_timestep, xlevs
    common /dynsoil_dble2D_1/ reg_thick,reg_x_surf,reg_tau_surf,reg_P_vol,reg_mean_age,reg_x_mean
    common /dynsoil_dble2D_2/ reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros
    common /dynsoil_dble2D_3/ reg_Li_Friv, reg_Li_Fsp, reg_Li_driv
    common /dynsoil_dble3D/   reg_z_prof,reg_tau_prof,reg_eros_galy_unit
    common /dynsoil_dble1D/   slope
    common /dynsoil_logi1D/   DS_missingpoints
    common /dynsoil_intg2D/   reg_ktop
    common /dynsoil_intvar/   DS_ofile_num
    common /dynsoil_charvar1/ DS_ofile_name
    common /dynsoil_charvar2/ DS_varout_name
    common /dynsoil_charvar3/ DS_varout_units, DS_varout_missvalname
    common /dynsoil_charvar4/ DS_restart_name
    common /dynsoil_dblevar/  DS_varout_missval
    ! vegetation:
    common /veget_int/ coupling_veget, icount_veget, ijump_veget, veget_occup
    common /veget/     veget_factor,veget_eros_factor
