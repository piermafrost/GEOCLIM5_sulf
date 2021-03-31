	implicit double precision (a-h,o-z)

	parameter (npixel=1920)
	parameter (nline=1200)
	parameter (lex=40)
	parameter (nidl_plot=7)

	dimension time(1:nline),T(1:nline,1:npixel),R(1:nline,1:npixel),area(1:npixel)
	dimension iline(1:nidl_plot),Tcont(1:nidl_plot,1:npixel)
	dimension Tsea(1:nidl_plot,1:npixel)

	open(1,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/T_map.200_eq1')
	open(2,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/R_map.200_eq1')
	open(3,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/var.200_eq1')
	open(4,file='/Users/salome/fortran/GEOCLIM2/pangea2014/Aire_200MaVegDyn.dat')
	
	open(100,file='T_map_continent.out')
	open(101,file='R_map.out')
	open(102,file='T_map_ocean.out')
	
	open(200,file='T_stage1_land.out')
	open(201,file='T_stage2_land.out')
	open(202,file='T_stage3_land.out')
	open(203,file='T_stage4_land.out')
	open(204,file='T_stage5_land.out')
	open(205,file='T_stage6_land.out')
	open(206,file='T_stage7_land.out')


	open(300,file='T_stage1_sea.out')
	open(301,file='T_stage2_sea.out')
	open(302,file='T_stage3_sea.out')
	open(303,file='T_stage4_sea.out')
	open(304,file='T_stage5_sea.out')
	open(305,file='T_stage6_sea.out')
	open(306,file='T_stage7_sea.out')

	do j=1,nline
	  read(1,*)(T(j,i),i=1,npixel)
	  read(2,*)(R(j,i),i=1,npixel)
	  read(3,*)time(j)
	end do

	do j=1,npixel
	  read(4,*)area(j)
	end do
	
	
	a=-1d+34
	do k=1,npixel
	  if (area(k).ne.0) then
	    write(100,*)T(lex,k)
	    write(101,*)R(lex,k)
	    write(102,*)a
	  else
	    write(100,*)a
	    write(101,*)a
	    write(102,*)T(lex,k)
	  endif
	end do


!  Preparing files for idl maps:
!  +++++++++++++++++++++++++++++	
	interval=int(nline/nidl_plot)

	do k=1,nidl_plot
	   iline(k)=0+k*interval
	   print*,iline(k)
	end do

	do j=1,nidl_plot
	 do k=1,npixel
	  if (area(k).ne.0) then
	    Tcont(j,k)=T(iline(j),k)
	    Tsea(j,k)=-1.d+34
	  else
	    Tcont(j,k)=-1.d+34
	    Tsea(j,k)=T(iline(j),k)
	  endif
	 end do
	end do
	
	
	write(200,13)(Tcont(1,k),k=1,npixel)
	write(201,13)(Tcont(2,k),k=1,npixel)
	write(202,13)(Tcont(3,k),k=1,npixel)
	write(203,13)(Tcont(4,k),k=1,npixel)
	write(204,13)(Tcont(5,k),k=1,npixel)
	write(205,13)(Tcont(6,k),k=1,npixel)
	write(206,13)(Tcont(7,k),k=1,npixel)

	write(300,13)(Tsea(1,k),k=1,npixel)
	write(301,13)(Tsea(2,k),k=1,npixel)
	write(302,13)(Tsea(3,k),k=1,npixel)
	write(303,13)(Tsea(4,k),k=1,npixel)
	write(304,13)(Tsea(5,k),k=1,npixel)
	write(305,13)(Tsea(6,k),k=1,npixel)
	write(306,13)(Tsea(7,k),k=1,npixel)
	
	

	
  13    format(1(1x,e15.5))

	
	stop
	end

