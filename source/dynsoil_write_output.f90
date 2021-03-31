module dynsoil_write_output_mod
implicit none

contains

subroutine dynsoil_write_output( filenum,filename,varname, nx,ny, t,temp,runoff,hsoil,xsurf,tausurf,z,tau, &
                                 reg_prod,reg_eros,reg_P_diss,reg_P_eros,xsurf_eros,x_mean,reg_mean_age,    &
                                 Li_Friv, Li_Fsp, Li_driv                                            )
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, &
                              put_var_real1D, put_var_real2D, put_var_real3D, put_var_real4D
  use netcdf
  integer, parameter:: nvar = 25
  integer, dimension(nvar), intent(in):: filenum
  character(len=*), dimension(nvar), intent(in):: filename, varname
  integer, intent(in):: nx, ny
  double precision, intent(in):: t
  double precision, dimension(:), intent(in):: temp,runoff
  double precision, dimension(:,:), intent(in):: hsoil,xsurf,tausurf,reg_prod,reg_eros,&
                                                 reg_P_diss,reg_P_eros,xsurf_eros,x_mean,reg_mean_age,&
                                                 Li_Friv,Li_Fsp,Li_driv
  double precision, dimension(:,:,:), intent(in):: z,tau
  integer, dimension(nvar):: nt, fileid, varid
  integer:: i,j,k,nlit,nlev,ierr,timevarid(1),dimid(1)

  !===================================== OUTPUT VARIABLES LIST: ======================================!
  !  OUTPUT: X, Y, litho, xlevs, t, area, litho_frac, slope, temp, runoff, h_soil, x_surf, tau_surf,  !
  !  var #:  1  2  3      4      5  6     7           8      9     10      11      12      13         !
  !          z, tau, Reg_prod, Reg_eros, reg_P_diss, reg_P_eros, x_surf_eros, x_P_mean, reg_mean_age, !
  !          14 15   16        17        18          19          20           21        22            !
  !          Li_Friv,  Li_Fsp,   Li_driv                                                              !
  !          23        24        25                                                                   !
  !===================================================================================================!


  ! get xlevels dimension length
  nlev = size(z,1)
  nlit = size(z,2)


!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!

    ! open file(s) and write time variable:
    !
    k = 0
    !
    do i = 7,8 ! open potential litho_frac and slope files
      if (filenum(i) > 0) then
        if ( filenum(i) > k ) then
          call open_file( filename(i) , j , mode=NF90_WRITE )
          fileid(i) = j
          call inquire_dim( fileid(i) , varname(5:5) , dimid )
          ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
          nt(i) = nt(i) + 1
          call inquire_var( fileid(i) , varname(5:5) , timevarid )
          call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
          k = filenum(i)
        else
          fileid(i) = fileid(filenum(i))
          nt(i) = nt(filenum(i))
        end if
      end if
    end do
    !
    do i = 9,nvar ! loop on all variables but litho_frac slope
      if (filenum(i) > 0) then
        if ( filenum(i) > k ) then
          call open_file( filename(i) , j , mode=NF90_WRITE )
          fileid(i) = j
          call inquire_dim( fileid(i) , varname(5:5) , dimid )
          ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
          nt(i) = nt(i) + 1
          call inquire_var( fileid(i) , varname(5:5) , timevarid )
          call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
          k = filenum(i)
        else
          fileid(i) = fileid(filenum(i))
          nt(i) = nt(filenum(i))
        end if
        call inquire_var( fileid(i) , varname(i:i) , varid(i:i) )
      end if
    end do


    !!!!!!!!!!!!!!!!!!!!
    ! write variables: !
    !!!!!!!!!!!!!!!!!!!!

    i = 9  ! TEMPERATURE
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) ,          real(reshape(temp, shape=(/nx,ny/))),  (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 10 ! RUNOFF
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) ,        real(reshape(runoff, shape=(/nx,ny/))),  (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 11 ! HSOIL
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,         real(reshape(hsoil, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 12 ! XSURF
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,         real(reshape(xsurf, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 13 ! TAUSURF
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,       real(reshape(tausurf, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 14 ! Z
    if (filenum(i)>0) then
      call put_var_real4D( fileid(i) , varid(i) ,             real(reshape(z, shape=(/nx,ny,nlev,nlit/), order=(/4,3,1,2/))) , &
                                                                                     (/1,1,1,1,nt(i)/) , (/nx,ny,nlit,nlev,1/) )
    end if
    i = 15 ! TAU
    if (filenum(i)>0) then
      call put_var_real4D( fileid(i) , varid(i) ,           real(reshape(tau, shape=(/nx,ny,nlev,nlit/), order=(/4,3,1,2/))) , &
                                                                                     (/1,1,1,1,nt(i)/) , (/nx,ny,nlit,nlev,1/) )
    end if
    i = 16 ! REG_PROD
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,      real(reshape(reg_prod, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 17 ! REG_EROS
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,     real(reshape(reg_eros, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 18 ! REG_P_DISS
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,    real(reshape(reg_P_diss, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 19 ! REG_P_EROS
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,    real(reshape(reg_P_eros, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 20 ! XSURF_EROS
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,    real(reshape(xsurf_eros, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 21 ! X_MEAN
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,        real(reshape(x_mean, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 22 ! REG_MEAN_AGE
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,  real(reshape(reg_mean_age, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 23 ! LI_FRIV
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,       real(reshape(Li_Friv, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 24 ! LI_FSP
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,        real(reshape(Li_Fsp, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 25 ! LI_DRIV
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) ,       real(reshape(Li_driv, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if

  !======================!
  ! Output file closing: !
  !======================!

  k = 0
  do i = 7,nvar
    if ( filenum(i) > k ) then
      call close_file(fileid(i))
      k = filenum(i)
    end if
  end do



end subroutine


end module
