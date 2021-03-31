'''
Export the netCDF data from ERA5 files to ascii files readable by GEOCLIM
(ravelled in 1D vector, with CO2 level before all values for temp and runf)
  Area: convert in 1e6km2
  Land area: convert in 1e6km2, multiply by land fraction ('lsm') and remove
             points with no slope or no lithology data (set land_area=0)
  Temperature: convert in deg Celsius
  Runoff: convert in cm/y
'''

import netCDF4 as nc
import numpy as np

PREINDUSTRIAL_CO2 = 286.7 #ppm

fa = nc.Dataset('grid_30minx30min_sphere.nc')
fl = nc.Dataset('landfrac_30minx30min.nc')
fc = nc.Dataset('ERA5_1981-2018_30minx30min.nc')
fs = nc.Dataset('../slope/slope_PD_SRTM_30min.nc')
flh = nc.Dataset('../lithology/lithology_fraction_6class_30minx30min.nc')

slopemask = fs['slope'][:,:].mask
lithomask = flh['lithfrac'][0,:,:].mask
area = 1e-12*fa['area'][:,:]
landarea = area*fl['lsm'][:,:]
landarea[np.logical_or(slopemask, lithomask)] = 0
temp = fc['t2m'][:,:] - 273.15
runoff = 100*fc['ro'][:,:]

np.savetxt('area_tot_1e6km2_ERA5.dat',  area.ravel(), fmt='%.7e')
np.savetxt('area_land_1e6km2_ERA5.dat', landarea.ravel(), fmt='%.7e')
np.savetxt('temperature_degC_ERA5.dat', np.concatenate(([286.7], temp.ravel())), fmt='%.7e')
np.savetxt('runoff_cmpy_ERA5.dat',      np.concatenate(([286.7], runoff.ravel())), fmt='%.7e')

fa.close()
fl.close()
fc.close()
fs.close()
flh.close()
