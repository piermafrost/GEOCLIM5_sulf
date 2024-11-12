import netCDF4 as nc

fT = nc.Dataset('/home/piermafrost/GitHub/GEOCLIM-DynSoil-steady-state/input/GFDL_temp_360_720_PD.nc')
fq = nc.Dataset('/home/piermafrost/GitHub/GEOCLIM-DynSoil-steady-state/input/GFDL_runoff_360_720_PD.nc')
fout = [nc.Dataset('GFDL_preind_qT_{:}ppm.nc'.format(int(c)), mode='w', format='NETCDF3_CLASSIC')
        for c in fT['lvl'][:].data]



for k,f in enumerate(fout):


    # Global Attribute
    f.title = 'GFDL CM2.0 preindustrial climatology, interpolated on a 0.5 degrees lon-lat grid.'
    f.CO2_level_ppm = fT['lvl'][k].data
    f.reference = 'T. L. Delworth et al., GFDL’s CM2 global coupled climate models. Part I: Formulation and simulation characteristics. J. Clim., 19, 643-674 (2006). DOI: 10.1175/JCLI3629.1    Park Y. et al., Emergence of the Southeast Asian islands as a driver for Neogene cooling. PNAS 117(41), 25319-25326 (2020). DOI: 10.1073/pnas.2011033117'

    # Copy dimensions and dim variables
    for dim in ['lon', 'lat']:
        f.createDimension(dim, fT.dimensions[dim].size)
        f.createVariable(dim, datatype=fT.variables[dim].datatype, dimensions=(dim,))
        # Attributes:
        for att in ['axis', 'units', 'long_name']:
            if hasattr(fT.variables[dim], att):
                setattr(f.variables[dim], att, getattr(fT.variables[dim], att))
        
        # Put variable
        f.variables[dim][:] = fT.variables[dim][:]

    # Copy main variables
    for fin,var in zip([fT,fq], ['tmp', 'rnf']):
        f.createVariable(var, datatype=fin.variables[var].datatype, dimensions=('lat', 'lon'), fill_value=fin.variables[var]._FillValue)
        # Attributes:
        for att in ['units', 'long_name']:
            if hasattr(fin.variables[var], att):
                setattr(f.variables[var], att, getattr(fin.variables[var], att))
        
        # Put variable
        f.variables[var][:,:] = fin.variables[var][k,:,:]

    f.close()

fT.close()
fq.close()

