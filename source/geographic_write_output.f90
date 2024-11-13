module geographic_write_output_mod
implicit none

contains

subroutine geographic_write_output(ofile_name, time_dimname, GEOG_outvar_info, t, cpvec, &
                                   temp, runoff, wth_allsil, wth_litho_wgh, wth_litho, ker_wth, bioC_exp, phos_wth)
  use netcdf
  use io_module,        only: netcdf_output_var
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, put_var
  include 'shape.inc' ! => nlon, nlat, nlitho
  include 'output_size.inc' ! => nGEOGoutvar
  character(len=*):: ofile_name, time_dimname
  type(netcdf_output_var), dimension(nGEOGoutvar), intent(in):: GEOG_outvar_info
  double precision, intent(in):: t
  double precision, dimension(5), intent(in):: cpvec
  double precision, dimension(:), intent(in):: temp,runoff, wth_allsil, ker_wth, bioC_exp, phos_wth 
  double precision, dimension(:,:), intent(in):: wth_litho_wgh, wth_litho
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

    do i = 3,7 ! => climatic parameters # 1 to 5
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real0D=real(cpvec(i-2)), stt=(/nt/), cnt=(/1/))
    end do
    !
    i = 8
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(temp, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 9
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(runoff, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 11
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(wth_allsil, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 12
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
 var_real3D=real(reshape(wth_litho_wgh, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 13
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
     var_real3D=real(reshape(wth_litho, shape=(/nlon,nlat,nlitho/), order=(/3,1,2/))), stt=(/1,1,1,nt/), cnt=(/nlon,nlat,nlitho,1/))
    !
    i = 14
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(ker_wth, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 15
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(bioC_exp, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))
    !
    i = 16
      if (GEOG_outvar_info(i)%writevar) &
        call put_var(fid, varname=GEOG_outvar_info(i)%vname, &
                     var_real2D=real(reshape(phos_wth, shape=(/nlon,nlat/))), stt=(/1,1,nt/), cnt=(/nlon,nlat,1/))


    ! Close output file
    ! -----------------

    call close_file(fid)



end subroutine


end module
