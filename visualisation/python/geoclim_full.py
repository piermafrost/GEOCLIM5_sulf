import numpy as np
import pylab as pl
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from matplotlib.backends.backend_pdf import PdfPages

global file_name,file_name2
global extension,path

# which simulation ?
path='/Users/yves/fortran/GEOCLIM4_thea/OUTPUT/'   # / at the end of the line
opath='/Users/yves/fortran/GEOCLIM4_thea/python/figs/'
extension=['200Ma_highCO2','200Ma_ref_highCO2']   #,'200Ma_peaks']     #,'090_evol2','090_evol3']



xaxis='time (Myr)'
nfile=len(extension)
file_name=list()     #empty list that will be filled with the names of files
file_name2=list()    #empty list that will be filled with the names of files
file_name3=list()    #empty list that will be filled with the names of files
file_name4=list()    #empty list that will be filled with the names of files
file_name5=list()    #empty list that will be filled with the names of files
file_name6=list()    #empty list that will be filled with the names of files
file_name7=list()    #empty list that will be filled with the names of files
file_name8=list()    #empty list that will be filled with the names of files
var_all=list()
amean=list()
asum=list()
fsize=tuple()
fsize=(17,14)


#-------------FUNCTIONS--------------------------
def making_file_name(j,i):
    index_box=['1','2','3','4','5','6','7','8','9']
    file_label='box'
    file_name=path+file_label+index_box[j-1]+'.'+extension[i-1]
    return file_name
#------------------------------------------------
def making_file_name2(i):
    file_label='var.'
    file_name2=path+file_label+extension[i-1]
    return file_name2
#------------------------------------------------
def making_file_name3(i):
    file_label='cflux.'
    file_name3=path+file_label+extension[i-1]
    return file_name3
#------------------------------------------------
def making_file_name4(j,i):
    index_box=['1','3','6','8']
    file_label='biodiv'
    file_name4=path+file_label+index_box[j-1]+'.'+extension[i-1]
    return file_name4
#------------------------------------------------
def making_file_name5(i):
    file_label='speciation_c13.'
    file_name5=path+file_label+extension[i-1]
    return file_name5
#------------------------------------------------
def making_file_name6(i):
    file_label='chimie.'
    file_name6=path+file_label+extension[i-1]
    return file_name6
#------------------------------------------------
def making_file_name7(i):
    file_label='lysocline.'
    file_name7=path+file_label+extension[i-1]
    return file_name7
#------------------------------------------------
def making_file_name8(i):
    file_label='species_count.'
    file_name8=path+file_label+extension[i-1]
    return file_name8
#------------------------------------------------






name=list()
namef=list()
namesp=list()
name=['DIC mol/m3','Alk eq/m3','PO4 mol/m3','Ca2+ mol/m3',
       'Sr2+ mol/m3','Sr PIC mol/m3','POP mol/m3','PIP mol/m3',
       'POC mol/m3','PIC mol/m3','O2 mol/m3','empty','DIC d13C',
       'PIC d13C','POC d13C','empty','seawater 87Sr/86Sr','PIC 87Sr/86Sr']
namef=['CO2 cons. granitic weath. (GtC/yr)','CO2 cons. basaltic weath. (GtC/yr)',
        'Carbonate weath. (GtC/yr)','Kerogen weath. (GtC/yr)','Reef deposition (GtC/yr)',
        'empty','Org C deposition (GtC/yr)']
namesp=['PNS species#','MLS species#','ES species#','PSS species#']
name_fi=['Cumul. CO2 cons granitic weath. (GtC)','Cumul. CO2 cons basaltic weath (GtC)',
         'Cumul. C released carbonate weath (GtC)','Cumul. C release ker weath (GtC)',
         'Cumul. C storage carb dep (GtC)','empty','Cumul. C storage org dep (GtC)']
ocean_basin=list()
ocean_basin_surface=list()
ocean_basin=['PNS','PND','MLS','MLT','MLD','ES','ED','PSS','PSD']
ocean_basin_surface=['PNS','MLS','ES','PSD']
styling=list()
styling=['-','--',':','-.','.','o','^','-','-']
ocean_volumes=tuple()
ocean_volumes=(48.36,130.14,28.7,258.3,772.3,0.8,0.98,27.65,74.4)



#Atmospheric gases------------------------------
print('atmospheric gases...')
plt.figure(num=1)
i=0
loop=True
while (loop):
  i=i+1
  file_name2.append( making_file_name2(i) )   #fill the list file_name with the name of the files
  var=np.loadtxt(file_name2[i-1])             #load files
  pl.subplot(2,2,3)
  plt.rcParams['xtick.labelsize'] = 9
  plt.rcParams['ytick.labelsize'] = 9
  pl.plot(var[:,0]/1e6,var[:,1],label=extension[i-1])             #plot O2
  pl.legend(loc='upper right',fontsize=9)
  pl.xlabel(xaxis,fontsize=9)
  pl.ylabel('atmospheric O2 (PAL)',fontsize=9)

  pl.subplot(2,1,1)
  plt.rcParams['xtick.labelsize'] = 9
  plt.rcParams['ytick.labelsize'] = 9
  pl.plot(var[:,0]/1e6,var[:,2]*280,label=extension[i-1])             #plot CO2
  pl.legend(loc='upper right',fontsize=7)
  pl.xlabel(xaxis,fontsize=9)
  pl.ylabel('atmospheric CO2 (ppm)',fontsize=9)
  
  Atm_flux=var[:,[4,5,6,7]]                                      #extracting the ocean-atm fluxes
  Atm_flux_total=list()                                          #calculating the sum of the ocean-atm fluxes
  Atm_flux_total=np.sum(Atm_flux,axis=1)                         #sum of the columns
  pl.subplot(2,2,4)
  plt.rcParams['xtick.labelsize'] = 7
  plt.rcParams['ytick.labelsize'] = 7
  pl.plot(var[:,0]/1e6,Atm_flux_total*12/(1e3*1e3*1e9),label=extension[i-1])
  pl.xlabel(xaxis,fontsize=9)
  pl.ylabel('total ocean-atm CO2 flux (GtC/yr)',fontsize=9)
  if  i>nfile-1:
     loop=False

#plt.show()
plt.tight_layout()
plt.savefig(opath+'atmospheric.pdf',format='pdf')
plt.close(1)



#Carbon fluxes------------------------------
print ('carbon fluxes...')
plt.figure(num=2,figsize=fsize)
i=0
loop=True
while (loop):
    j=0
    i=i+1
    file_name3.append( making_file_name3(i) )   #fill the list file_name with the name of the files
    var=np.loadtxt(file_name3[i-1])             #load files
    loop2=True

    while (loop2):
      j=j+1
      pl.subplot(3,3,j)
      plt.rcParams['xtick.labelsize'] = 12
      plt.rcParams['ytick.labelsize'] = 12
      pl.plot(var[:,0]/1e6,var[:,j]*12/(1e3*1e3*1e9),label=extension[i-1])             #plot global fluxes
      pl.legend(loc='upper right',fontsize=12)
      pl.xlabel(xaxis,fontsize=12)
      pl.ylabel(namef[j-1],fontsize=12)
      if j>6:
          loop2=False
    if i>nfile-1:
      loop=False

#plt.show()
plt.tight_layout()
plt.savefig(opath+'Carbon_fluxes.pdf',format='pdf')
plt.close(2)






plt.figure(num=3,figsize=fsize)
print('calculating integrated fluxes...')
i=0
loop=True
while (loop):
    loop2=True
    i=i+1
    file_name3.append( making_file_name3(i) )   #fill the list file_name with the name of the files
    var=np.loadtxt(file_name3[i-1])             #load files
    integrated_flux=np.zeros((len(var[:,0]),7))
    timespan=var[1,0]-var[0,0]
    
    kplot=0    
    while (loop2):
      kplot=kplot+1
      for k in range(len(var[:,0])):
        j=0
        loop_sum=True
        while (loop_sum):
          integrated_flux[k,kplot-1]=integrated_flux[k,kplot-1]+timespan*(var[j,kplot]-var[0,kplot])*12/(1e3*1e3*1e9)
          j=j+1
          if j > k:
            loop_sum=False

      pl.subplot(3,3,kplot)            
      plt.rcParams['xtick.labelsize'] = 12
      plt.rcParams['ytick.labelsize'] = 12
      pl.plot(var[:,0]/1e6,integrated_flux[:,kplot-1],label=extension[i-1])
      pl.legend(loc='upper right',fontsize=12)
      pl.xlabel(xaxis,fontsize=12)
      pl.ylabel(name_fi[kplot-1],fontsize=12)
      if kplot>6:
          loop2=False      

    if i>nfile-1:
      loop=False

print ('...done')
    
#    pl.subplot(3,3,9)
#    plt.rcParams['xtick.labelsize'] = 9
#    plt.rcParams['ytick.labelsize'] = 9
#    pl.plot(var[:,0]/1e6,var[:,22]/1e12,label=extension[i-1])             #plot global fluxes
#    pl.plot(var[:,0]/1e6,var[:,23]/1e12,'-',label=extension[i-1])             #plot global fluxes
#    pl.legend(loc='upper right',fontsize=9)
#    pl.xlabel(xaxis,fontsize=8)
#    pl.ylabel('ES organic carbon burial (Tmol/yr)',fontsize=9)

#    pl.subplot(3,3,9)
#    plt.rcParams['xtick.labelsize'] = 9
#    plt.rcParams['ytick.labelsize'] = 9
#    pl.plot(var[:,0]/1e6,(var[:,18]+var[:,21]+var[:,25])/1e12,label=extension[i-1])             #plot global fluxes
#    pl.legend(loc='upper right',fontsize=9)
#    pl.xlabel(xaxis,fontsize=8)
#    pl.ylabel('Open ocean organic carbon burial (Tmol/yr)',fontsize=9)

#plt.show()
plt.tight_layout()
plt.savefig(opath+'Integrated_fluxes.pdf',format='pdf')
plt.close(3)








##Biomass plots-------------------------------------
#plt.figure(num=3,figsize=fsize)
##i = file index
##j = basin index
##l = variable index
#i=0
#loop=True
#while (loop):
#    j=0
#    i=i+1
#    loop2=True
#    while (loop2):
#        j=j+1
#        #        print k,i,j
#        l=0
#        file_name4.append( making_file_name4(j,i) )   #fill the list file_name with the name of the files
#        k=(i-1)*4+j-1                              #list index (this is a vector, not a matrix
#        #        print k,file_name4[k]
#        biodiv=np.loadtxt(file_name4[k])               #load files
#        nlength=len(biodiv[0,:])
#        biodiv_total=np.sum(biodiv,axis=1)             #sum of the columns
#        pl.subplot(2,2,j)
#        plt.rcParams['xtick.labelsize'] = 12
#        plt.rcParams['ytick.labelsize'] = 12
#        pl.plot(var[:,0]/1e6,biodiv_total,label=extension[i-1])  #,linestyle=styling[j-1])
#        pl.legend(loc='upper right',fontsize=13)
#        pl.xlabel(xaxis,fontsize=12)
#        pl.ylabel('Total biomass (mol C) '+ocean_basin_surface[j-1],fontsize=12)
#        pl.yscale('log')
#        if j>3:
#           loop2=False
#    if i>nfile-1:
#      loop=False
#
#plt.savefig(opath+'biodiv.pdf',format='pdf')
#plt.close(3)
#
#
#
###Species number plots-------------------------------------
#plt.figure(num=7,figsize=fsize)
#i=0
#loop=True
#while (loop):
#    loop2=True
#    i=i+1
#    file_name8.append( making_file_name8(i))   #fill the list file_name with the name of the files
#    varsp=np.loadtxt(file_name8[i-1])             #load files
#    j=0
#    while (loop2):
#        j=j+1
#        pl.subplot(2,2,j)
#        pl.bar(varsp[:,0],varsp[:,j], facecolor='#9999ff', edgecolor='white',label=extension[i-1])
#        pl.plot(varsp[:,0],varsp[:,j],label=extension[i-1])             #plot number of species
#        pl.ylim(0,100)
#        pl.legend(loc='upper right',fontsize=9)
#        pl.xlabel(xaxis)
#        pl.ylabel(namesp[j-1])
#        plt.rcParams['xtick.labelsize'] = 10
#        plt.rcParams['ytick.labelsize'] = 10
#        if j>3:
#             loop2=False
#        
#    if  i>nfile-1:
#      loop=False
#
##plt.show()
#plt.savefig(opath+'species.pdf',format='pdf')
#plt.close(7)
#


##Biomass diversity mapsplots-------------------------------------
#plt.figure(num=4,figsize=fsize)
##i = file index
##j = basin index
##l = variable index
#i=0
#loop=True
#while (loop):
#    j=0
#    i=i+1
#    loop2=True
#    while (loop2):
#        j=j+1
#        print k,i,j
#        l=0
#        file_name4.append( making_file_name4(j,i) )   #fill the list file_name with the name of the files
#        k=(i-1)*4+j-1                              #list index (this is a vector, not a matrix
#        print k,file_name4[k]
#        biodiv=np.loadtxt(file_name4[k])               #load files
#        nlength=len(biodiv[0,:])
#        pl.subplot(2,2,j)
#        biomass_div=np.log10(biodiv)
#        pl.contourf (biomass_div,cmap='jet')
#        plt.rcParams['xtick.labelsize'] = 12
#        plt.rcParams['ytick.labelsize'] = 12
#        pl.xlabel('species #',fontsize=12)
#        pl.ylabel('Time (Myr)',fontsize=12)
#        if j>3:
#            loop2=False
#    if i>nfile-1:
#        loop=False
#
#plt.savefig('biodiv_map.pdf',format='pdf')
#plt.close(4)




#Ocean chemistry------------------------------
print ('ocean chemistry...')
plt.figure(num=5,figsize=fsize)
i=0
loop=True
while (loop):
        i=i+1
        file_name6.append( making_file_name6(i))   #fill the list file_name with the name of the files
        chemical=np.loadtxt(file_name6[i-1])             #load files
        pl.subplot(2,2,1)
        pl.plot(chemical[:,0]/1e6,chemical[:,10],label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('seawater ES pH')
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10
        
        pl.subplot(2,2,2)
        pl.plot(chemical[:,0]/1e6,(chemical[:,7]-chemical[:,9]),label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('seawater MLS-MLD pH gradient')
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10
        
        pl.subplot(2,2,3)
        pl.plot(chemical[:,0]/1e6,chemical[:,19],label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('seawater ES omega')
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10

        pl.subplot(4,2,6)
        pl.plot(chemical[:,0]/1e6,chemical[:,16],label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('seawater MLS omega')
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10

        pl.subplot(4,2,8)
        pl.plot(chemical[:,0]/1e6,chemical[:,18],label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('seawater MLD omega')
        plt.rcParams['xtick.labelsize'] = 10
        plt.rcParams['ytick.labelsize'] = 10

        if  i>nfile-1:
            loop=False

#plt.show()
plt.savefig(opath+'chemical.pdf',format='pdf')
plt.close(5)



#lysocline------------------------------
print ('lysocline...')
plt.figure(num=6,figsize=fsize)
i=0
loop=True
while (loop):
        i=i+1
        file_name7.append( making_file_name7(i))   #fill the list file_name with the name of the files
        lysocline=np.loadtxt(file_name7[i-1])             #load files
        pl.subplot(3,1,1)
        pl.plot(lysocline[:,0]/1e6,-lysocline[:,5]*1e3,label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('MLD lysocline depth (m)')
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13

        pl.subplot(3,1,2)
        pl.plot(lysocline[:,0]/1e6,-lysocline[:,2]*1e3,label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('PND lysocline depth (m)')
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13

        pl.subplot(3,1,3)
        pl.plot(lysocline[:,0]/1e6,-lysocline[:,9]*1e3,label=extension[i-1])             #plot pH
        pl.legend(loc='upper right',fontsize=9)
        pl.xlabel(xaxis)
        pl.ylabel('PNS lysocline depth (m)')
        plt.rcParams['xtick.labelsize'] = 13
        plt.rcParams['ytick.labelsize'] = 13


        if  i>nfile-1:
            loop=False

#plt.show()
plt.savefig(opath+'lysocline.pdf',format='pdf')
plt.close(6)




#Box plots-------------------------------------
print ('oceanic reservoirs...')
pdf_pages=PdfPages(opath+'boxes.pdf')

store_c13=np.zeros((len(lysocline),21,nfile))

#plt.figure(num=2,figsize=fsize)
#i = file index
#j = basin index
#l = variable index
i=0
loop=True
while (loop):
    fig = plt.figure(figsize=fsize)
    j=0
    i=i+1
    loop2=True
    while (loop2):
        j=j+1
        l=0
        file_name.append( making_file_name(j,i) )   #fill the list file_name with the name of the files
        k=(i-1)*9+j-1                              #list index (this is a vector, not a matrix
        box=np.loadtxt(file_name[k])               #load files
        storage=(i-1)*9+5
        print(store_c13.shape)
        print(k,file_name[k],box.shape)
        print('storage',storage)
        if k == storage:
            print('storing')
            store_c13[:,:,i-1]=box[:,:]                          #storing the box6 for use in 13C plots
        loop3=True
        while (loop3):
            l=l+1
            plt.rcParams['xtick.labelsize'] = 6
            plt.rcParams['ytick.labelsize'] = 6
            pl.subplot(5,5,l)
            pl.plot(box[:,0]/1e6,box[:,l],label=ocean_basin[j-1])  #,linestyle=styling[j-1])
            pl.legend(loc='upper right',fontsize=6)
            pl.xlabel(xaxis,fontsize=8)
            pl.ylabel(name[l-1],fontsize=8)
            pl.title(extension[i-1])
            if l>17:
               loop3=False
        if j>8:
            loop2=False
    plt.tight_layout()
    pdf_pages.savefig(fig)
    if i>nfile-1:
      loop=False

pdf_pages.close()
#plt.savefig('boxes.pdf',format='pdf')
#plt.close(2)

#plt.show()



#Seawater d13C------------------------------
print ('carbon isotopes...')
plt.figure(num=4,figsize=fsize)
i=0
loop=True
while (loop):
        i=i+1
        file_name5.append(making_file_name5(i))   #fill the list file_name with the name of the files
        var_dc=np.loadtxt(file_name5[i-1])             #load files
        plt.rcParams['xtick.labelsize'] = 12
        plt.rcParams['ytick.labelsize'] = 12
        pl.subplot(2,2,1)
        pl.plot(var_dc[:,0]/1e6,var_dc[:,6]*1000,label='H$_2$CO$_3$ '+extension[i-1])             #plot D13C
        pl.plot(var_dc[:,0]/1e6,var_dc[:,15]*1000,label='HCO$_3$ '+extension[i-1])             #plot D13C
        pl.plot(var_dc[:,0]/1e6,var_dc[:,24]*1000,label='CO$_3$ '+extension[i-1])             #plot D13C
        pl.legend(loc='upper right',fontsize=12)
        pl.xlabel(xaxis,fontsize=15)
        pl.ylabel('epicontinental surface $\delta^{13}$C (permil)',fontsize=15)

        pl.subplot(2,2,2)
        pl.plot(var_dc[:,0]/1e6,(store_c13[:,13,i-1])*1000,label='DIC '+extension[i-1])             #plot D13C
        pl.plot(var_dc[:,0]/1e6,(store_c13[:,14,i-1])*1000,label='PIC '+extension[i-1])             #plot D13C
        pl.legend(loc='upper right',fontsize=12)
        pl.xlabel(xaxis,fontsize=15)
        pl.ylabel('epicontinental surface $\delta^{13}$C (permil)',fontsize=15)

        pl.ylabel('epicontinental surface $\delta^{13}$C (permil)',fontsize=15)
        pl.subplot(2,2,3)
        pl.plot(var_dc[:,0]/1e6,(store_c13[:,15,i-1])*1000,label='POC '+extension[i-1])             #plot D13C
        pl.legend(loc='upper right',fontsize=12)
        pl.xlabel(xaxis,fontsize=15)
        pl.ylabel('epicontinental surface $\delta^{13}$C (permil)',fontsize=15)

        if  i>nfile-1:
            loop=False

#plt.show()
plt.savefig(opath+'d13C.pdf',format='pdf')
plt.close(4)



