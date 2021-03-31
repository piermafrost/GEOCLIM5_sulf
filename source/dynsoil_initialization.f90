module dynsoil_initialization_mod
implicit none

contains


subroutine dynsoil_initialization( xlevs, slope,area, h,xsurf,tausurf,z,tau, reg_prod,reg_eros,reg_P_diss,reg_P_eros,xsurf_eros,  &
                                   x_mean,reg_mean_age,reg_P_vol, ktop, Li_Friv,Li_Fsp,Li_driv, slopemissval,missingpoints )
  use dynsoil, only: get_ktop_and_Pvol
  use io_module, only: check_landcells, check_landcells_generic, check_invalid

  !===================================== OUTPUT VARIABLES LIST: ======================================!
  !  OUTPUT: X, Y, litho, xlevs, t, area, litho_frac, slope, temp, runoff, h_soil, x_surf, tau_surf,  !
  !  var #:  1  2  3      4      5  6     7           8      9     10      11      12      13         !
  !          z, tau, Reg_prod, Reg_eros, reg_P_diss, reg_P_eros, x_surf_eros, x_P_mean, reg_mean_age, !
  !          14 15   16        17        18          19          20           21        22            !
  !          Li_Friv,  Li_Fsp,   Li_driv                                                              !
  !          23        24        25                                                                   !
  !===================================================================================================!

  double precision, intent(in), dimension(:):: xlevs
  double precision, intent(inout), dimension(:):: area, slope
  double precision, intent(inout), dimension(:,:):: h, xsurf, tausurf
  double precision, intent(inout), dimension(:,:,:):: z, tau
  double precision, intent(out), dimension(:,:):: reg_prod, reg_eros, reg_P_diss, reg_P_eros, xsurf_eros, &
                                                  x_mean, reg_mean_age, reg_P_vol, Li_Friv, Li_Fsp, Li_driv
  integer, intent(out), dimension(:,:):: ktop
  double precision, intent(in):: slopemissval
  logical, intent(inout), dimension(:):: missingpoints

  integer:: i, npxl

  ! Global GEOCLIM variable:
  integer, dimension(5):: ERROR_HANDLING_OPTION
  common /error/ ERROR_HANDLING_OPTION
  ! ERROR_HANDLING_OPTION tells what to do for missing points issue (2) and invalid slope (3)


  !============================= PRELIMINARY =============================!

  ! get size of data:
  npxl = size(area)


  ! Check consistency between missing-points in GEOCLIM inputs (ie: continental area = 0)
  ! in dynsoil inputs (missingpoint = .true.) and in slope

  call check_landcells('slope', slopemissval, area, ERROR_HANDLING_OPTION(2), var1D=slope)
  call check_landcells_generic('WARNING: found missing value on continental cells of DynSoil init. variables (h, x_s or tau_s)', &
                               missingpoints, area, ERROR_HANDLING_OPTION(2))

  where (area==0) missingpoints=.true.

  ! Check for null slope
  call check_invalid('slope', area, ERROR_HANDLING_OPTION(3), var1D=slope)


  !============================= COMPUTATION =============================!

  ! get ktop and cation mass on each points:
  call get_ktop_and_Pvol( xlevs ,xsurf,h,z, missingpoints, ktop,reg_P_vol )



end subroutine



end module
