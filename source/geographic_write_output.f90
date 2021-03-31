module geographic_write_output_mod
implicit none

contains

subroutine geographic_write_output( filenum, filename, varname,  nx, ny, t, &
                                    temp, runoff, wth_allsil, wth_litho_wgh, wth_litho, ker_wth, bioC_exp, phos_wth )
  use netcdf_io_module, only: open_file, close_file, inquire_var, inquire_dim, put_var_real1D, put_var_real2D, put_var_real3D
  use netcdf
  integer, parameter:: nvar = 15
  integer, dimension(nvar), intent(in):: filenum
  character(len=*), dimension(nvar), intent(in):: filename, varname
  integer, intent(in):: nx, ny
  integer:: nlit
  double precision, intent(in):: t
  double precision, dimension(:), intent(in):: temp,runoff, wth_allsil, ker_wth, bioC_exp, phos_wth 
  double precision, dimension(:,:), intent(in):: wth_litho_wgh, wth_litho
  integer, dimension(nvar):: nt, fileid, varid
  integer:: i,j,k,ierr,timevarid(1),dimid(1)

  !===================================== OUTPUT VARIABLES LIST: ===================================================================!
  !  OUTPUT: lon, lat, litho time, area, litho_frac temperature, runoff, slope, sil_weathering, weighted_weathering                !
  !  var #:  1    2    3     4     5     6          7            8       9      10              11                                 !
  !          unweighted_weathering, ker_weathering, bio_Corg_export, phos_weathering                                               !
  !          12                     13              14               15                                                            !
  !================================================================================================================================!



    nlit = size(wth_litho,1)

!=================================================================================================!
!------------------------------------- OPENNING AND WRITING: -------------------------------------!
!=================================================================================================!

    ! open file(s) and write time variable:
    !
    k = 0
    !
    i = 6 ! open potential litho_frac file
    if (filenum(i) > 0) then
      if ( filenum(i) > k ) then
        call open_file( filename(i) , j , mode=NF90_WRITE )
        fileid(i) = j
        call inquire_dim( fileid(i) , varname(4:4) , dimid )
        ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
        nt(i) = nt(i) + 1
        call inquire_var( fileid(i) , varname(4:4) , timevarid )
        call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
        k = filenum(i)
      else
        fileid(i) = fileid(filenum(i))
        nt(i) = nt(filenum(i))
      end if
    end if
    !
    do i = 7,8 ! loop on all first variables
      if (filenum(i) > 0) then
        if ( filenum(i) > k ) then
          call open_file( filename(i) , j , mode=NF90_WRITE )
          fileid(i) = j
          call inquire_dim( fileid(i) , varname(4:4) , dimid )
          ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
          nt(i) = nt(i) + 1
          call inquire_var( fileid(i) , varname(4:4) , timevarid )
          call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
          k = filenum(i)
        else
          fileid(i) = fileid(filenum(i))
          nt(i) = nt(filenum(i))
        end if
        call inquire_var( fileid(i) , varname(i:i) , varid(i:i) )
      end if
    end do
    !
    i = 9 ! open potential slope file
    if (filenum(i) > 0) then
      if ( filenum(i) > k ) then
        call open_file( filename(i) , j , mode=NF90_WRITE )
        fileid(i) = j
        call inquire_dim( fileid(i) , varname(4:4) , dimid )
        ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
        nt(i) = nt(i) + 1
        call inquire_var( fileid(i) , varname(4:4) , timevarid )
        call put_var_real1D( fileid(i), timevarid(1), (/real(t)/) , nt(i:i) , (/1/) )
        k = filenum(i)
      else
        fileid(i) = fileid(filenum(i))
        nt(i) = nt(filenum(i))
      end if
    end if
    !
    do i = 10,nvar ! loop on all variables after slope
      if (filenum(i) > 0) then
        if ( filenum(i) > k ) then
          call open_file( filename(i) , j , mode=NF90_WRITE )
          fileid(i) = j
          call inquire_dim( fileid(i) , varname(4:4) , dimid )
          ierr = nf90_inquire_dimension( fileid(i) , dimid(1) , len=nt(i) )
          nt(i) = nt(i) + 1
          call inquire_var( fileid(i) , varname(4:4) , timevarid )
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

    i = 7
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(temp, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 8
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(runoff, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 10
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(wth_allsil, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 11
    if (filenum(i)>0) then
      call put_var_real3D( fileid(i) , varid(i) , real(reshape(wth_litho_wgh, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                     (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 12
    if (filenum(i)>0) then
    call put_var_real3D( fileid(i) , varid(i) , real(reshape(wth_litho, shape=(/nx,ny,nlit/), order=(/3,1,2/))) , &
                                                                                (/1,1,1,nt(i)/) , (/nx,ny,nlit,1/) )
    end if
    i = 13
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(ker_wth, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 14
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(bioC_exp, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
    end if
    i = 15
    if (filenum(i)>0) then
      call put_var_real2D( fileid(i) , varid(i) , real(reshape(phos_wth, shape=(/nx,ny/))) , (/1,1,nt(i)/) , (/nx,ny,1/) )
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
