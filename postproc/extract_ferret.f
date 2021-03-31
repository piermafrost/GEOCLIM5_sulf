	implicit double precision (a-h,o-z)

	parameter (npixel=1920)
	parameter (nline=4000)
	parameter (lex=40)
	parameter (nselect=2)

	dimension time(1:nline),T(1:nline,1:npixel),R(1:nline,1:npixel),area(1:npixel)
	dimension Tcont(1:npixel),Rcont(1:npixel)
	dimension Tsea(1:npixel)

	open(1,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/T_map.200_deg')
	open(2,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/R_map.200_deg')
	open(3,file='/Users/salome/fortran/GEOCLIM2/OUTPUT/var.200_deg')
	open(4,file='/Users/salome/fortran/GEOCLIM2/pangea2014/Aire_200MaVegDyn.dat')
	
	open(100,file='T_map_continent2.out')
	open(101,file='R_map2.out')
	open(102,file='T_map_ocean2.out')
	
	
	do j=1,nline
	read (1,*) (T(j,k),k=1,npixel)
	read (2,*) (R(j,k),k=1,npixel)
	end do
	
	do k=1,npixel
	read (4,*) area(k)
	end do
		
	do k=1,npixel
	if (area(k).eq.0) then 
	Tcont(k)=-1d+34
	Tsea(k)=T(nselect,k)
	Rcont(k)=-1d+34
	else 
	Tcont(k)=T(nselect,k)
	Tsea(k)=-1d+34
	Rcont(k)=R(nselect,k)
	end if
	end do
	
	do k=1,npixel
	write (100,*) Tcont(k)
	write (102,*) Tsea(k)
	write (101,*) Rcont(k)
	end do
	
	stop
	end
	  