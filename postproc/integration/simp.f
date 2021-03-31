      implicit double precision(a-h,o-z)
      parameter(dim1=2000)
	  parameter(Fod_ss=1.60501e+12)

      common /trans/ cfl(dim1),time(dim1)
	  common /integrating/ cfl_ready(dim1)
      common /trans2/ tsample(1000)
	  

      open(unit=5,file='/Users/yves/fortran/GEOCLIM4_thea/OUTPUT/cflux.200Ma_peaks')

      open(unit=8,file='integr.res')

      do j=1,dim1
          read(5,*)time(j),not,not,not,not,not,not,cfl(j)
		  cfl_ready(j)=(cfl(j)-Fod_ss)*12/(1e3*1e3*1e9)
      end do
	  

c	  do j=1,801
c	    call interd(tim(j),conc_int,tp,conc,20,k1,k2)
c	    write(8,*)tim(j),xmet(j),conc_int   !discharge in m3/s
c	  end do

	  sint=0.
	  t0=0.
	  ni=1000000
	  
	  do j=1,1000
		  tsample(j)=1000.*j
	  end do
	  
	  do j=1,1000
            call simps(t0,tsample(j),sint,ni)
            write(8,*)tsample(j),sint
	end do

      stop
      end
      

      SUBROUTINE SIMPS (X0,XN,SINT,NI)
C     --------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	  
      LI=2*NI
      H=(XN-X0)/FLOAT(LI)
      SUM=0.
      DO 100 I=0,LI
       K=2*INT((FLOAT(I)+0.01)/2.)
       IF(K.NE.I)THEN
           FAC=4.
       ELSEIF((I.EQ.0).OR.(I.EQ.LI))THEN
           FAC=1.
         ELSE
           FAC=2.
       ENDIF
       X=X0+FLOAT(I)*H
       CALL FCT(X,FX)
       SUM=SUM+FAC*FX
 100  CONTINUE
      SINT=SUM*H/3.
      RETURN
      END

      subroutine fct (t,fx)
c     ---------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter(dim1=2000)
      common /integrating/ cfl_ready(dim1)
      common /trans2/ tsample(1000)
      common /trans/ cfl(dim1),time(dim1)

   	  call interd(t,res,time,cfl_ready,dim1,k1,k2)

      fx=res
c	  print*,res,k1,k2
	  
	  
      return
      end



c	+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine interd(t,x,tdat,xdat,nn,k1,k2)
c     -----------------------------------------
      implicit double precision (a-h,o-z)
      dimension tdat(nn),xdat(nn)

 
      do 200 k=1,nn
       if(t.le.tdat(k))then
         k2=k
         go to 210
       endif
 200  continue
      k2=nn
 210  continue
      k1=k2-1
      x=xdat(k1)+(xdat(k2)-xdat(k1))*(t-tdat(k1))/
     1         (tdat(k2)-tdat(k1))
      return
      end

    
