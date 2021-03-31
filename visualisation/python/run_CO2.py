import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages

global file_name
global extension,path

path='/Users/yves/fortran/GEOCLIM4/calibration/PD/'   # / at the end of the line
opath='/Users/yves/fortran/GEOCLIM4/python/figs/'
extension=('out')     #,'090_evol2','090_evol3']
xaxis='time (Myr)'





styling=list()
styling=['-','--',':','-.','.','o','^','-','-']



#runoff------------------------------
plt.figure(num=1)
file_label='run'
var=np.loadtxt(path+file_label+'.'+extension)             #load files

pl.contourf(var)

#pl.subplot(2,2,3)
#i=0
#loop=True
#while (loop):
#    i=i+1
#    pl.plot(var[:,0]/1e6,var[:,5000+i])
#    if i>5000:
#       loop=False



#pl.legend(loc='upper right',fontsize=9)
#pl.xlabel(xaxis)
#pl.ylabel('atmospheric O2 (PAL)')
#plt.rcParams['xtick.labelsize'] = 7
#plt.rcParams['ytick.labelsize'] = 7
#pl.legend(loc='upper right',fontsize=7)
#pl.xlabel(xaxis)
#pl.ylabel('atmospheric CO2 (PAL)')
  
#runoff=var[:,[1:16384]]                                      #extracting the ocean-atm fluxes
#print type(runoff)
#  Atm_flux_total=list()                                          #calculating the sum of the ocean-atm fluxes
#  Atm_flux_total=np.sum(Atm_flux,axis=1)                         #sum of the columns
#  pl.subplot(2,2,4)
#  plt.rcParams['xtick.labelsize'] = 7
#  plt.rcParams['ytick.labelsize'] = 7
#  pl.plot(var[:,0]/1e6,Atm_flux_total/1e12,label=extension[i-1])
#  pl.xlabel(xaxis)
#  pl.ylabel('total ocean-atm CO2 flux (Tmol/yr)',fontsize=9)
#  if  i>nfile-1:
#     loop=False

plt.show()
#plt.savefig(opath+'atmospheric.pdf',format='pdf')
plt.close(1)





