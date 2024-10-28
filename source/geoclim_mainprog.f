!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   + COupled Model of BIogeochemical cycles aNd climatE +
!   + COMBINE                                            +
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++

!   Calibrated with ERA5 reanalysis (continental temperature and runoff),
!   STRM (slope) and Hartmann & Moosdorf 2012 (lithology map),
!   resolution: 0.5deg x 0.5deg


program geoclim



    use omp_lib
    use constante,                     only: PI_n_CO2_atm
    use GCM_io_module,                 only: load_climatology
    use geoclim_create_output_mod,     only: geoclim_create_output
    use geographic_create_output_mod,  only: geographic_create_output
    use dynsoil_read_input_mod,        only: dynsoil_read_input
    use dynsoil_create_output_mod,     only: dynsoil_create_output
    use dynsoil_initialization_mod,    only: dynsoil_initialization
    use io_module,                     only: check_invalid
    use output_netcdf2ascii_mod,       only: output_netcdf2ascii
    use utils,                         only: read_comment, set_error_handling_option
    implicit none
    include 'combine_foam.inc'
    integer:: ierr





!   ####################################################################################
    open(unit=0, status='old', action='read', file=geoclim_path//'config/IO_CONDITIONS')
!   ####################################################################################


    ! Set the variable ERROR_HANDLING_OPTION telling the code what to do when
    ! some errors are encountered
    call set_error_handling_option()


    ! read run name, output directory and open cond_p20:
    call READ_MAIN_IO_CONDITIONS(0)


    ! Display run name
    print *
    print *
    n = len(trim(run_name))
    do k = 1,n+13
        write(*, fmt='(A1)', advance='no') '#'
    end do
    print *
    write(*, fmt='(A)') '##  RUN: '//trim(run_name)//'  ##'
    do k = 1,n+13
        write(*, fmt='(A1)', advance='no') '#'
    end do
    print *


    call read_conditions
    tbegin=tsta
    tend=tfi
    stept=ts
!    check=xjump
    time=tbegin
    compteur=int(xjump)+1
    htry=ts
    eps=1.d-8
    icount_geogprint = ijump_geogprint
    icount_veget = ijump_veget
    if (coupling_veget==0) then
        !no modification
        veget_factor = 1
        veget_eros_factor = 1
    end if


!   ############ FILE OPENNING ############ !
    call OPEN_ASCII_FILES(0)
!   ####################################### !


!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   Read continental inputs:

    if (land_input_mode == 'ascii') then

        ref_x_axis = 0d0
        ref_y_axis = 0d0
        ! => tells the code to create them later

        do i=1,npixel  !S ---> N
            read(7,*)areaclimber(i)  !continental area
        end do
        close(unit=7)

        do i=1,npixel  !S ---> N
            read(301,*)areaEarth(i)  !continental and oceanic area
        end do
        close(unit=301)

        ! Note: both areas are expected to be in m2

        do j=1,nclimber
            read(30,*)co2climber(nclimber+1-j), &
                    (Tairclimber(i,nclimber+1-j),i=1,npixel)
            read(31,*)dummy,(Runclimber(i,nclimber+1-j),i=1,npixel)
        end do
        close(unit=30)
        close(unit=31)

        GMSTclimber = -1d99 ! no definition of GMST in inputs


    elseif (land_input_mode == 'GCM') then

        call load_climatology(333)
        ! Get global variables:
        !    'co2climber'
        !    'ref_x_axis'
        !    'ref_y_axis'
        !    'areaEarth'
        !    'areaclimber'
        !    'Tairclimber'
        !    'Runclimber'
        !    'GMSTclimber'

    end if

    ! Units conversion:
    !
    ! area: m2 -> 1e6km2:
    areaclimber = areaclimber/1e12
    areaEarth = areaEarth/1e12
    !
    ! CO2: ppmv -> PAL
    co2climber = co2climber/280.

    ! get total land area
    areatot = sum(areaclimber)

    call load_lithology(304)
    ! get global variable 'litho_frac'

!   Check for negative runoff
    print *
    call check_invalid('runoff', areaclimber, ERROR_HANDLING_OPTION(3), var2D=Runclimber, axis=2)


!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        


    ! Read oceanic inputs

    call read_oceanic_temperature(32, co2climber, Toceclimber)

    ! Store GMST in last box (atmosphere) of Toceclimber:
    Toceclimber(nbasin,:) = GMSTclimber

    ! Obsolete: hypso file
    !do j=1,13
    !    read(4,*)deephyp(j),shyp(j)
    !end do
    !close(unit=4)

    areal=0.

    call basin_geometry()

    close(unit=333)

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   <><><><><><><><><><><><>
!   read initial conditions
!   <><><><><><><><><><><><>


!   COMBINE:
!   --------

    do j=1,nvar
        read(2,*)y(j)
    enddo
    close(unit=2)

    ! extensive variables: concentration => chemical amount
    do j = 1,11*nbasin
        y(j) = y(j)*vol(j)
    end do
    do j = 1+19*nbasin,nvar
        y(j) = y(j)*vol(j)
    end do

    if (nclimber==1) then
        ! fixed atm CO2 run => set atmospheric CO2 as defined by the unique CO2 level
        ! = co2_level (PAL) * Pre-Ind moles of CO2 in atm (0.0508d18)
        ! (atm PCO2 is 12th variable of last basin (nbasin))
        y(12*nbasin) = co2climber(1)*PI_n_CO2_atm
    endif

    call varset(y,nvar)


!   ECOLOGICAL NETWORK MODULE (=> LEFT FOR FUTURE DEVELOPMENT):
!   -----------------------------------------------------------

    if (coupling_ecogeo) then
        print *, 'ECOGEO module not available'
        stop -1
    end if


!    DynSoil
!    -------

    if (coupling_dynsoil) then
        ! Note: there are 12 uncommented lines to read:
        !  - name of restart file to create
        !  - initialization mode
        !  - initialization file name
        !  - 8 initialization variables
        !  - slope file and variable name

        call read_comment(0)
        read(unit=0,fmt=*) dummychar, DS_restart_name
        call read_comment(0)
        read(unit=0,fmt=*) DS_init_mode
        ! read init file info or create init cond
        call dynsoil_read_input( DS_init_mode, 0,                                                 &
                                 ref_x_axis,ref_y_axis,xlevs,                                     &
                                 reg_thick,reg_x_surf,reg_tau_surf,reg_z_prof,reg_tau_prof,slope, &
                                 DS_missingpoints, DS_slopemissval                                )

        ! initialisation of rest of dynsoil variables:
        call dynsoil_initialization( xlevs, slope,areaclimber,                                  &
                                     reg_thick,reg_x_surf,reg_tau_surf,reg_z_prof,reg_tau_prof, &
                                     reg_prod,reg_eros,reg_P_diss,reg_P_eros,reg_x_surf_eros,   &
                                     reg_x_mean,reg_mean_age,reg_P_vol, reg_ktop,               &
                                     reg_Li_Friv,reg_Li_Fsp,reg_Li_driv,                        &
                                     DS_slopemissval,DS_missingpoints                           )

        DS_timestep = ijump_DS_integration * ijump_cont_weath * ts
        icount_DS_int = ijump_DS_integration
        icount_DS_pri = ijump_DS_print

    else

        ! skip the 12 uncommented lines
        do i = 1,12
            call read_comment(0)
            read(unit=0,fmt=*)
        end do

    end if

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   Create geographic coordinate (longitude and latitude) if they were not given in input data
    if (all(ref_x_axis==0d0) .and. all(ref_y_axis==0d0)) then
        ref_x_axis = (/ ( ((dble(i)-0.5)/nlon) * 360 - 180 + 360/nlon/2   , i=1,nlon ) /)
        ref_y_axis = (/ ( ((dble(j)-0.5)/nlat) * 180 - 90  + 180/nlat/2   , j=1,nlat ) /)
    end if

!   Get the list of continental pixels (areaclimber > 0)
    call get_cont_pixel()




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   <><><><><><><><><><>
!   Create output files
!   <><><><><><><><><><>

!   read combine (oceanic) output conditions and create output file(s):
    call geoclim_create_output( 0, output_path,run_name,vol,oce_surf,surf_sedi,                 &
                                GEO_ofile_name,GEO_ofile_num,GEO_varout_name,GEO_varout_missval )
    ! Special case of GMST (optional input)
    where(GMSTclimber==-1d99)
        GMSTclimber           = GEO_varout_missval(95)
        Toceclimber(nbasin,:) = GEO_varout_missval(95)
    end where

!   read geographic output conditions and create output file(s):
    call geographic_create_output( 0, output_path,run_name,ref_x_axis,ref_y_axis,areaclimber,slope, &
                 litho_frac,GEOG_ofile_name,GEOG_ofile_num,GEOG_varout_name,GEOG_varout_missval     )

!   read dynsoil output conditions and create output file(s):
    if (coupling_dynsoil) then
        call dynsoil_create_output( 0, output_path,run_name,                                    &
                                    xlevs,ref_x_axis,ref_y_axis,areaclimber,litho_frac,slope,   &
                                    DS_ofile_num,DS_ofile_name,DS_varout_name,DS_varout_units,  &
                                    DS_varout_missvalname,DS_varout_missval                     )

    else
        do i = 1,nDSvar
            call read_comment(0)
            read(0,*)
        end do
    end if

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        


!   Restart file not yet created:
    fog = 0


!   close IO condition file:
    close(0)

!   signal that the current run is done with the 2 configuration files,
!   they can be safely modified.
    open(unit=0, file=geoclim_path//'.config-queue', status='old', action='readwrite', iostat=ierr)
    if (ierr/=0)  open(unit=0, file=geoclim_path//'.config-queue', status='replace', action='readwrite')
    ierr=0
    do while (ierr==0) ! read until end of file
        read(unit=0, fmt=*, iostat=ierr)
    end do
    if (ierr>0) read(unit=0, fmt=*) ! if error other than end-of-file raised, re-execute reading action to display it
    backspace(unit=0)
    write(unit=0, fmt='(A)') 'free (run: '//trim(run_name)//')'
    close(unit=0)




!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   Put missing-value on non-continental points (area==0) of time-varying geographic output variables
!   (fillvalues are read during creation of output files, time-constant variables already written)

    do j = 1,npixel
        if (areaclimber(j)==0) then
            Tclim(j)           = GEOG_varout_missval(7)
            Runclim(j)         = GEOG_varout_missval(8)
            wth_allsil(j)      = GEOG_varout_missval(10)
            wth_litho_wgh(:,j) = GEOG_varout_missval(11)
            wth_litho(:,j)     = GEOG_varout_missval(12)
            fker(j)            = GEOG_varout_missval(13)
            POC_export_rate(j) = GEOG_varout_missval(14)
            fp(j)              = GEOG_varout_missval(15)
        end if
    end do

    ! DynSoil output variables
    if (coupling_dynsoil) then
        do j = 1,npixel
            ! missing-value for initial fluxes on ALL pixels (can only be computed between two time steps)
            reg_prod(:,j)        = DS_varout_missval(16)
            reg_eros(:,j)        = DS_varout_missval(17)
            reg_P_diss(:,j)      = DS_varout_missval(18)
            reg_P_eros(:,j)      = DS_varout_missval(19)
            reg_x_surf_eros(:,j) = DS_varout_missval(20)
            ! missing-value on non-continental pixels for other output variables
            if (areaclimber(j)==0) then
                reg_thick(:,j)      = DS_varout_missval(11)
                reg_x_surf(:,j)     = DS_varout_missval(12)
                reg_tau_surf(:,j)   = DS_varout_missval(13)
                reg_z_prof(:,:,j)   = DS_varout_missval(14)
                reg_tau_prof(:,:,j) = DS_varout_missval(15)
                reg_x_mean(:,j)     = DS_varout_missval(21)
                reg_mean_age(:,j)   = DS_varout_missval(22)
                reg_Li_Friv(:,j)    = DS_varout_missval(23)
                reg_Li_Fsp(:,j)     = DS_varout_missval(24)
                reg_Li_driv(:,j)    = DS_varout_missval(25)
            else
                ! missing-value above top of regolith for 3D variables
                do k = 1,nlitho
                    reg_z_prof(   reg_ktop(k,j)+1:nlitho, k, j ) = DS_varout_missval(14)
                    reg_tau_prof( reg_ktop(k,j)+1:nlitho, k, j ) = DS_varout_missval(15)
                end do
            end if
        end do
    end if

!   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        

!   End of initialization:

    htot=0.
    icompteur=int(xjump)+1
    icount_cont_weath = ijump_cont_weath

!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!   ++++++++++++++++++++++++++++++++++++++++++++++++++++++        




!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   !
!   dummy integration (do not actualized
!   prognostic variables) to get initiale fluxes

    print *
    print *
    print *, 'Compute initial values of fluxes (for outputs only)'

    ! SAVE INITIAL VAR:
    y_0 = y
    ! dynsoil prognostic variables:
    if (coupling_dynsoil) then
      reg_thick_0    = reg_thick
      reg_x_surf_0   = reg_x_surf
      reg_tau_surf_0 = reg_tau_surf
      reg_P_vol_0    = reg_P_vol
      reg_z_prof_0   = reg_z_prof
      reg_tau_prof_0 = reg_tau_prof
      reg_ktop_0     = reg_ktop
    end if

    call read_veget
    call creades(time)
    call rk4(y,dydx,nvar,time,htry,yout)

    do j=1,nvar
      y(j)=yout(j)
    end do


    ! RESET VARIABLES WITH INITIAL VALUES:
    y = y_0
    if (coupling_dynsoil) then
      reg_thick    = reg_thick_0
      reg_x_surf   = reg_x_surf_0
      reg_tau_surf = reg_tau_surf_0
      reg_P_vol    = reg_P_vol_0
      reg_z_prof   = reg_z_prof_0
      reg_tau_prof = reg_tau_prof_0
      reg_ktop     = reg_ktop_0
      icount_DS_int = ijump_DS_integration
    end if
    call varset(y,nvar)
    icount_cont_weath = ijump_cont_weath
    icount_veget = ijump_veget

    call printf(time,icompteur,y)

    ! For not printing the initial state a second time but still counting it:
    if (coupling_dynsoil) icount_DS_pri = -1
    icount_geogprint = -1
    compteur=-1
    icompteur=-1

!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   !




!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!   %%         TIME LOOP         %%
!   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    print *
    print *
    print *, 'Start run'
    call date_and_time(values=computer_time)
    write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
        computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
        computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)


!   &&&&&&&&&&  start iteration  &&&&&&&&&&
    do while (time.le.tend)

      call read_veget()
      call creades(time)
      call rk4(y,dydx,nvar,time,htry,yout)
      do j=1,nvar
        y(j)=yout(j)
      end do

      call varset(y,nvar)

      call printf(time,icompteur,y)
      time=time+htry

    end do
!   &&&&&&&&&&  end iteration  &&&&&&&&&&


    ! call 'printf' again only for creating restart files (in case asked to do at end of run)
    icompteur = -1
    icount_geogprint = -1
    icount_DS_pri = -1
    ! => do not print anything
    ageYprint = time - 2*ts
    ! => create restart, if not already created
    call printf(time, icompteur, y)


    print *
    print *, 'End of run'
    call date_and_time(values=computer_time)
    write(*, fmt='(I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I3.3)') &
        computer_time(1),'/',computer_time(2),'/',computer_time(3),' ', &
        computer_time(5),':',computer_time(6),':',computer_time(7),':',computer_time(8)


    ! CONVERT NETCDF OUTPUTS TO ASCII OUTPUTS:
    if (convert2ascii==1) call output_netcdf2ascii()





end program

