module dynsoil_initialization_mod
implicit none

contains

subroutine dynsoil_initialization( xlevs, h, xsurf, z, missingpoints, reg_P_vol, ktop )
  use dynsoil, only: get_ktop_and_Pvol

  double precision, intent(in), dimension(:):: xlevs
  double precision, intent(in), dimension(:,:):: h, xsurf
  double precision, intent(in), dimension(:,:,:):: z
  double precision, intent(out), dimension(:,:):: reg_P_vol
  integer, intent(out), dimension(:,:):: ktop
  logical, intent(in), dimension(:):: missingpoints

  ! get ktop and cation mass on each points:
  call get_ktop_and_Pvol( xlevs ,xsurf,h,z, missingpoints, ktop,reg_P_vol )

end subroutine

end module
