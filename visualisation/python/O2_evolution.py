import netCDF4 as nc
from matplotlib import pyplot as plt

f = nc.Dataset('../../OUTPUT/geoclim_output.test.nc')

# CO2 and O2 time evolution
plt.plot(f['time'], f['CO2_atm_level'])
plt.plot(f['time'], f['O2_atm_level'])
plt.legend(['CO2', 'O2'])
plt.xlabel('time (yr)')
plt.ylabel('atmospheric level (PAL)')
plt.show()

# OC burial time evolution
plt.plot(f['time'], f['org_C_tot_dep_flux'])
plt.xlabel('time (yr)')
plt.ylabel('OC burial (mol/yr)')
plt.show()
