import netCDF4 as nc
import numpy as np

# Orbital parameters time-series
# ------------------------------

ecc = lambda t: 0.006 + (0.05-0.006)*(1-np.cos(2*np.pi*t/405e3))/2
obl = lambda t: 22.1 + (24.5-22.1)*(1-np.cos(2*np.pi*t/41e3))/2
pre = lambda t: 360*(t/23e3)

# Field GEOCLIM inputs (climate and slope)
# ----------------------------------------

# Default fillvalue:
FILLVAL = -1e36

LON = np.array([0., 1., 2.])
LAT = np.array([-0.5, 0.5])
SHP = (LAT.size, LON.size) 
MASK = np.ones(SHP, dtype=bool)
# Define continental points:
MASK[0,1] = False
MASK[1,2] = False

AREA = (510e12/np.prod(SHP))*np.ones(SHP)
LAND = np.ones(SHP, dtype=AREA.dtype)
LAND[MASK] = 0.

SLOPE = 0.0123*np.ones(SHP, dtype='float32')
SLOPE[MASK] = FILLVAL
SLOPE = np.ma.masked_where(MASK, SLOPE)

# Default units
UNITS = {'area': 'm2', 'landfrac': '-', 'temperature': 'degrees_celsius', 'runoff': 'cm/yr', 'slope': 'm/m'}

def netcdf_output(fname, **kwargs):
    '''
    function to save fields (keyword arguments) in the netCDF file fname
    '''

    fout = nc.Dataset(fname, mode='w', data_model='NETCDF3_CLASSIC')

    fout.createDimension('lon', size=LON.size)
    fout.createDimension('lat', size=LAT.size)

    fout.createVariable('lon', datatype=LON.dtype, dimensions=('lon',))
    fout.variables['lon'].setncattr('units', '-')
    fout.variables['lon'][:] = LON
    fout.createVariable('lat', datatype=LAT.dtype, dimensions=('lat',))
    fout.variables['lat'].setncattr('units', '-')
    fout.variables['lat'][:] = LAT
    for var in kwargs:
        fout.createVariable(var, datatype=kwargs[var].dtype, dimensions=('lat','lon'), fill_value=FILLVAL)
        fout.variables[var].setncattr('units', UNITS[var])
        fout.variables[var][:] = kwargs[var]

    fout.close()


# Create GEOCLIM input fields for all orbital configuration
# =========================================================

Temp = lambda x: x*np.ones(SHP, dtype='float32')
Runf = lambda x: np.ma.masked_where(MASK, x*np.ones(SHP, dtype='float32'))

netcdf_output('grid.nc', area=AREA, landfrac=LAND)
netcdf_output('slope.nc', slope=SLOPE)
netcdf_output('clim_ecc-low_obl-low_pre-0.nc',   temperature=Temp(13),   runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-low_pre-90.nc',  temperature=Temp(13),   runoff=Runf(32.4))
netcdf_output('clim_ecc-low_obl-low_pre-180.nc', temperature=Temp(13),   runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-low_pre-270.nc', temperature=Temp(13),   runoff=Runf(32.6))
netcdf_output('clim_ecc-hig_obl-low_pre-0.nc',   temperature=Temp(15),   runoff=Runf(32.5))
netcdf_output('clim_ecc-hig_obl-low_pre-90.nc',  temperature=Temp(15),   runoff=Runf(30.))
netcdf_output('clim_ecc-hig_obl-low_pre-180.nc', temperature=Temp(15),   runoff=Runf(32.7))
netcdf_output('clim_ecc-hig_obl-low_pre-270.nc', temperature=Temp(15),   runoff=Runf(35.))
netcdf_output('clim_ecc-low_obl-mid_pre-0.nc',   temperature=Temp(12.8), runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-mid_pre-90.nc',  temperature=Temp(12.8), runoff=Runf(32.4))
netcdf_output('clim_ecc-low_obl-mid_pre-180.nc', temperature=Temp(12.8), runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-mid_pre-270.nc', temperature=Temp(12.8), runoff=Runf(32.6))
netcdf_output('clim_ecc-hig_obl-mid_pre-0.nc',   temperature=Temp(14.8), runoff=Runf(32.5))
netcdf_output('clim_ecc-hig_obl-mid_pre-90.nc',  temperature=Temp(14.8), runoff=Runf(30.))
netcdf_output('clim_ecc-hig_obl-mid_pre-180.nc', temperature=Temp(14.8), runoff=Runf(32.7))
netcdf_output('clim_ecc-hig_obl-mid_pre-270.nc', temperature=Temp(14.8), runoff=Runf(35.))
netcdf_output('clim_ecc-low_obl-hig_pre-0.nc',   temperature=Temp(12.6), runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-hig_pre-90.nc',  temperature=Temp(12.6), runoff=Runf(32.4))
netcdf_output('clim_ecc-low_obl-hig_pre-180.nc', temperature=Temp(12.6), runoff=Runf(32.5))
netcdf_output('clim_ecc-low_obl-hig_pre-270.nc', temperature=Temp(12.6), runoff=Runf(32.6))
netcdf_output('clim_ecc-hig_obl-hig_pre-0.nc',   temperature=Temp(14.6), runoff=Runf(32.5))
netcdf_output('clim_ecc-hig_obl-hig_pre-90.nc',  temperature=Temp(14.6), runoff=Runf(30.))
netcdf_output('clim_ecc-hig_obl-hig_pre-180.nc', temperature=Temp(14.6), runoff=Runf(32.7))
netcdf_output('clim_ecc-hig_obl-hig_pre-270.nc', temperature=Temp(14.6), runoff=Runf(35.))


# Create input file of orbital parameter time-series
# ==================================================

TEND = 850e3
DT = 1e3

t = np.arange(0., TEND, DT)
np.savetxt('time.dat', t, fmt='%.7e', delimiter='\n')
np.savetxt('orbit_params.dat', np.array([obl(t), ecc(t), pre(t)]).transpose(),
           fmt='%.7f', delimiter=', ', header='# obliquity, eccentricity, precession')

