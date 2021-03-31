pro interp

close,1,2,3

set_plot,'PS'
DEVICE,FILE='/Users/yves/fortran/GEOCLIM3/jnl/interp.ps'


read_mat,'/Users/yves/fortran/GEOCLIM3/calibration/265Ma/Area265MaCO.dat',area

ngrid=16384
nco2=5

null_val=0.

temp=fltarr(ngrid+1,nco2)
run=fltarr(ngrid+1,nco2)
openr,1,'/Users/yves/fortran/GEOCLIM3/calibration/265Ma/temp.out'
openr,2,'/Users/yves/fortran/GEOCLIM3/calibration/265Ma/run.out'

readf,1,temp
readf,2,run
close,1,2

ktot=indgen(1)
k=indgen(1)
CO2=fltarr(1,nco2)



for j=0,nco2-1 do begin
    CO2(j)=temp(0,j)
end

ktot=0
for j=0,ngrid-1 do begin
;    if (area(j) ne 0) then begin
       ktot=ktot+1
;    end
end
print,'continental grid cell number',ktot

temp_plot=fltarr(ktot,nco2)
run_plot=fltarr(ktot,nco2)
A=fltarr(ktot)
Ar=fltarr(ktot)
B=fltarr(ktot)
Br=fltarr(ktot)
Aparam=fltarr(1)
Bparam=fltarr(1)
Arparam=fltarr(1)
Brparam=fltarr(1)

k=0

for j=0,ngrid-1 do begin
;    if (area(j) ne 0) then begin
	  k=k+1
	  for i=0,nco2-1 do begin
	   temp_plot(k-1,i)=temp(j+1,i)
	   run_plot(k-1,i)=run(j+1,i)
	  end
;	end
end


!P.multi=[0,3,3]
!X.thick=3
!Y.thick=3
plot,co2(0,*)/280,temp_plot(1,*),psym=3,xtitle='atmospheric CO!i2!n (PAL)',$
     ytitle='air temperature FOAM',yrange=[min(temp_plot),max(temp_plot)],ystyle=1
for j=1,ktot-1 do begin
    oplot,co2(0,*)/280,temp_plot(j,*),psym=3
end

plot,temp_plot(1,*),run_plot(1,*),psym=3,xtitle='air temperature FOAM',$
     ytitle='runoff FOAM',xrange=[min(temp_plot),max(temp_plot)],xstyle=1,$
     yrange=[min(run_plot),150],ystyle=1
for j=1,ktot-1 do begin
    oplot,temp_plot(j,*),run(j,*),psym=3
end





for j=0,ktot-1 do begin
    result=linfit(alog(CO2/280),temp_plot(j,*))
    A(j)=result(0)
    B(j)=result(1)
    result=linfit(temp_plot(j,*),run_plot(j,*))
    Ar(j)=result(0)
    Br(j)=result(1)
end


;=========================================================================
;plotting the family of parameters fot the T vs CO2 correlation
plot,A,B,psym=3,/ystyle,xtitle='parameter A',ytitle='parameter B',$
                 title='A+B log(CO!i2!n)'


result=linfit(A,B)
Aparam=result(0)
Bparam=result(1)

oplot,A,Aparam+Bparam*A,thick=2
;=========================================================================






;=========================================================================
;plotting the family of parameters fot the runoff vs T correlation
plot,Ar,Br,psym=3,/ystyle,xtitle='parameter Ar',ytitle='parameter Br',$
                 title='Ar+Br T'


result=linfit(Ar,Br)
Arparam=result(0)
Brparam=result(1)

oplot,Ar,Arparam+Brparam*Ar,thick=2
;=========================================================================



;=========================================================================
pix=100
plot,CO2/280,A(pix)+B(pix)*alog(CO2/280),xtitle='atmospheric CO!i2!n (PAL)',$
             ytitle='air temperature (C)',title='grid element test'   ;direct log fit
oplot,CO2/280,temp_plot(pix,*),psym=1,symsize=2  
oplot,CO2/280,A(pix)+(Aparam+Bparam*A(pix))*alog(CO2/280),psym=6
x=findgen(80)+1
oplot,x,A(pix)+(Aparam+Bparam*A(pix))*alog(x),line=1

pix=3000
plot,CO2/280,A(pix)+B(pix)*alog(CO2/280),xtitle='atmospheric CO!i2!n (PAL)',$
             ytitle='air temperature (C)',title='grid element test'   ;direct log fit
oplot,CO2/280,temp_plot(pix,*),psym=1,symsize=2  
oplot,CO2/280,A(pix)+(Aparam+Bparam*A(pix))*alog(CO2/280),psym=6
x=findgen(80)+1
oplot,x,A(pix)+(Aparam+Bparam*A(pix))*alog(x),line=1
;=========================================================================



;=========================================================================
pix=100
plot,temp_plot(pix,*),Ar(pix)+Br(pix)*temp_plot(pix,*),xtitle='air temperature (C)',$
             ytitle='runoff (cm/yr)',title='grid element test'   ;direct fit
oplot,temp_plot(pix,*),run_plot(pix,*),psym=1,symsize=2  
oplot,temp_plot(pix,*),Ar(pix)+(Arparam+Brparam*Ar(pix))*temp_plot(pix,*),psym=6
x=findgen(60)-20
oplot,x,Ar(pix)+(Arparam+Brparam*Ar(pix))*x,line=1

pix=3000
plot,temp_plot(pix,*),Ar(pix)+Br(pix)*temp_plot(pix,*),xtitle='air temperature (C)',$
             ytitle='runoff (cm/yr)',title='grid element test'   ;direct fit
oplot,temp_plot(pix,*),run_plot(pix,*),psym=1,symsize=2  
oplot,temp_plot(pix,*),Ar(pix)+(Arparam+Brparam*Ar(pix))*temp_plot(pix,*),psym=6
x=findgen(60)-20
oplot,x,Ar(pix)+(Arparam+Brparam*Ar(pix))*x,line=1
;=========================================================================





print,'parameters:'
print,'+++++++++++'
print,'A= ',Aparam
print,'B= ',Bparam


openw,3,'/Users/yves/fortran/GEOCLIM3/jnl/factorT.out'
increment=0
for j=0,ngrid-1 do begin
;    if (area(j) eq 0) then begin
;        printf,3,null_val,null_val
;    endif else begin
        printf,3,A(increment),B(increment)
        increment=increment+1
;    endelse
end

close,3
openw,3,'/Users/yves/fortran/GEOCLIM3/jnl/factorR.out'
increment=0
for j=0,ngrid-1 do begin
;    if (area(j) eq 0) then begin
;        printf,3,null_val,null_val
;    endif else begin
        printf,3,Ar(increment),Br(increment)
        increment=increment+1
;    endelse
end

close,3

!P.multi=[0,0,0]

device,/close

end
