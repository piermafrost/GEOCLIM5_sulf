module dynsoil_write_output_mod
implicit none

contains

subroutine dynsoil_write_output(ofile_name, time_dimname, DYNS_outvar_info, t, temp, runoff, hsoil, xsurf, tausurf, z, tau,        &
                             reg_prod, reg_eros, reg_P_diss, reg_P_eros, xsurf_eros, x_mean, reg_mean_age, Li_Friv, Li_Fsp, Li_driv)
                                 
  use io_module, only: netcdf_output_var
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, put_var
  use netcdf
  include 'shape.inc' ! => nlon, nlat, nlitho, nDSlev
  include 'output_size.inc' ! => nDYNSoutvar
  character(len=*):: ofile_name, time_dimname
  type(netcdf_output_var), dimension(nDYNSoutvar), intent(in):: DYNS_outvar_info
  double precision, intent(in):: t
  double precision, dimension(:), intent(in):: temp, runoff
  double precision, dimension(:,:), intent(in):: hsoil, xsurf, tausurf, reg_prod, reg_eros, &
                                                 reg_P_diss, reg_P_eros, xsurf_eros, x_mean, reg_mean_age, &
                                                 Li_Friv, Li_Fsp, Li_driv
  double precision, dimension(:,:,:), intent(in):: z, tau
  integer:: i, ierr, nt, fid, timevarid, dimid


!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!


    ! Open file, get size of time dimension, and put current time
    ! -----------------------------------------------------------
    call open_file(ofile_name, fid, mode=NF90_WRITE)
    call inquire_dim(fid, time_dimname, dimid)
    ierr = nf90_inquire_dimension(fid, dimid, len=nt)
    nt = nt + 1
    call inquire_var(fid, time_dimname, timevarid)
    call put_var(fid, timevarid, var_real0D=real(t), stt=(/nt/), cnt=(/1/))



    !<><><><><><><><><><><><><><><><><>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<> %  write output variables  % <>!
    !<> %%%%%%%%%%%%%%%%%%%%%%%%%%%% <>!
    !<><><><><><><><><><><><><><><><><>!

    ! Note: the following list of blocks needs to be updated if one wants to add new output variables
    ! ***********************************************************************************************

    i = 4  ! TEMPERATURE
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real2D=real(reshape(temp, shape=(/nlon,nlat/))), &
                   stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 5 ! RUNOFF
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real2D=real(reshape(runoff, shape=(/nlon,nlat/))), &
                   stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 6 ! HSOIL
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(hsoil, shape=(/nlon,nlat,nlitho/), &
                                                       order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 7 ! XSURF
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(xsurf, shape=(/nlon,nlat,nlitho/), &
                                                       order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 8 ! TAUSURF
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(tausurf, shape=(/nlon,nlat,nlitho/), &
                                                         order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 9 ! Z
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real4D=real(reshape(z, shape=(/nlon,nlat,nDSlev,nlitho/), &
                                               order=(/4,3,1,2/))), stt=(/1,1,1,1,nt/), cnt=(/nlon,nlat,nlitho,nDSlev,1/))
    !
    i = 10 ! TAU
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real4D=real(reshape(tau, shape=(/nlon,nlat,nDSlev,nlitho/), &
                                                 order=(/4,3,1,2/))), stt=(/1,1,1,1,nt/), cnt=(/nlon,nlat,nlitho,nDSlev,1/))
    !
    i = 11 ! REG_PROD
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(reg_prod, shape=(/nlon,nlat,nlitho/), &
                                                          order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 12 ! REG_EROS
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(reg_eros, shape=(/nlon,nlat,nlitho/), &
                                                          order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 13 ! REG_P_DISS
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(reg_P_diss, shape=(/nlon,nlat,nlitho/), &
                                                            order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 14 ! REG_P_EROS
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(reg_P_eros, shape=(/nlon,nlat,nlitho/), &
                                                            order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 15 ! XSURF_EROS
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(xsurf_eros, shape=(/nlon,nlat,nlitho/), &
                                                            order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 16 ! X_MEAN
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(x_mean, shape=(/nlon,nlat,nlitho/), &
                                                        order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 17 ! REG_MEAN_AGE
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(reg_mean_age, shape=(/nlon,nlat,nlitho/), &
                                                              order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 18 ! LI_FRIV
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(Li_Friv, shape=(/nlon,nlat,nlitho/), &
                                                         order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 19 ! LI_FSP
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(Li_Fsp, shape=(/nlon,nlat,nlitho/), &
                                                        order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 20 ! LI_DRIV
    if (DYNS_outvar_info(i)%writevar) &
      call put_var(fid, varname=DYNS_outvar_info(i)%vname, var_real3D=real(reshape(Li_driv, shape=(/nlon,nlat,nlitho/), &
                                                         order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))


    ! Close output file
    ! -----------------

    call close_file(fid)



end subroutine


end module
