## Example of use:
##
## write_oceanic_temp(co2=[560, 840, 1120, 1210, 1305, 1400, 1680, 1960, 2240, 4480, 8960],
##                    input_files=['560ppm.nc', '840ppm.nc', '1120ppm.nc', '1210ppm.nc', '1305ppm.nc',
##                                 '1400ppm.nc', '1680ppm.nc', '1960ppm.nc', '2240ppm.nc', '4480ppm.nc', '8960ppm.nc'],
##                    root='/home/piermafrost/Downloads/150Ma/merg150_VegDef_AdjCSol_EccN_ocean_',
##                    latitude='lat', z='lev', temperature='TEMP',
##                    y_weight='lw', z_weight='thickness')

import netCDF4 as nc
import numpy as np
from units import temperature_units, length_units, latitude_units




###############################
# SPECIFIC GEOCLIM PARAMETERS #
###############################

# Oceanic basin definition
NBASIN = 9 # not counting atmosphere box
H_EPICONT = 200  # (m) maximum ocean floor depth of the epicontinental box
H_SURF    = 100  # (m) depth of surface boxes (ie, starting depth of thermocline or epicont deep box)
H_THERMO  = 1000 # (m) depth of thermocline boxes (ie, starting depth of deep box)
MIDLAT_RANGE = [-60, 60] # (degrees N): latitudinal range of inner ocean (ie, excluding high-latitude oceans)

# Arbitrary parameters
ATM_VOL = 1.
ATM_SURF = 0.363e15
ATM_SEDI_SURF = 0.357e15

# Technical parameters (MUST BE KEPT UNMODIFIED)
SURF_CONVERSION_FACTOR = 1e-12 # surface expressed in 1e9km2
VOLUME_CONVERSION_FACTOR = 1e-15 # volume expressed in 1e6km3


##########################
# MASK OF GEOCLIM BASINS #
##########################

def geoclim_basin_mask(nav_lat, nav_z, nav_depth):
    '''
    Create a mask telling if points are in (T) or outside (F) each of the 9 GEOCLIM
    basins, for a given latitude array "nav_lat", local depth array "nav_z",
    and ocean floor depth array "nav_depth".
    the returned mask will rank-4, the 1st dimension being GEOCLIM basins, the last
    3 being the dimensions of the 3-D input fields (i.e., latitude, z and depth)

    nav_lat, nav_z and nav_depth MUST BE RANK-3, even if they are not defined on
    all dimensions (e.g., latitude and floor depth should be defined on horizontal,
    dimensions, whereas z should be defined on the vertical dimension). Use
    degenerated (size-1) dimensions for those "extra" dimensions.
    The 3 dimensions must obviously correspond (ie, be in the same order) to the
    dimensions of 3D oceanic variables (e.g., temperature).

    The GEOCLIM basin definition must be the following:
      0. N high-lat, surface
      1. N high-lat, deep (incl. thermo)
      2. inner oce, surface
      3. inner oce, thermocline
      4. inner oce, deep
      5. epicont, surface
      6. epicont, deep
      7. S high-lat, surface
      8. S high-lat, deep (incl. thermo)
    '''


    shp = np.maximum(np.maximum(nav_lat.shape, nav_z.shape), nav_depth.shape)

    mask = np.ones(np.concatenate(([NBASIN], shp)), dtype=bool)


    epicont_mask = (nav_depth <= H_EPICONT)
    north_mask   = (nav_lat > MIDLAT_RANGE[1])
    south_mask   = (nav_lat < MIDLAT_RANGE[0])

    # Northern high-lat basins
    for i in [0,1]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], north_mask)

    # Inner ocean basins
    for i in [2,3,4]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], np.logical_and(~north_mask, ~south_mask))

    # Southern high-lat basins
    for i in [7,8]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], south_mask)

    # Epicontinental basins
    for i in [5,6]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], epicont_mask)

    # Non epicontinental basins
    for i in [0,1,2,3,4,7,8]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], ~epicont_mask)

    # Surface basins
    for i in [0,2,5,7]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], nav_z<=H_SURF)

    # Thermocline, high-lat deep and epicontinental deep basins
    for i in [1,3,6,8]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], nav_z>H_SURF)

    # Thermocline basin
    for i in [3]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], nav_z<=H_THERMO)

    # Inner ocean deep basin
    for i in [4]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], nav_z>H_THERMO)


    return mask


#######################
# AUXILIARY FUNCTIONS #
#######################

def find_depth(var_mask, z, dim=None):
    '''
    Return a 2D (horizontal) array giving the seafloor depth of each point.
    The seafloor depth is determined as follows:
      1. Find the index the last 'False' of 3D boolean array "var_mask" in the vertical
         (z) dimension (yield a 2D "horizontal" index array `idx`).
      2. Pick the value of "z" (1D "vertical" array) for each "bottom index"
         (i.e., return `z[idx]`). If all the value of var_mask in one column are 'True',
         the depth of that point will be '0'
    "var_mask" must be a rank-3 boolean array. "z" must be a rank-1 array, or a rank-3
    array with 2 degenerated (size-1) dimensions.
    use optional argument "dim" to specified which of "var_mask" dimensions is the
    vertical one.
    If the vertical (z) dimension is not specified, the algorithm will:
      1- consider the only non-degenerated dimension of "z", if it is rank-3
      2- consider the only dimension of "var_mask" that have the same length
         than "z"
    If the procedure yield a unique consistent possibility, the user will be notified,
    of the automatic decision. If not (for instance, if several dimensions of "var_mask"
    have the same length than "z", an error will be raised).
    '''

    if var_mask.ndim != 3:
        raise ValueError('"var_mask" must be rank-3')

    # find vertical dimension
    if z.ndim == 3:
        if dim is None:
            dim = np.argwhere(np.array(z.shape)>1)
            if dim.size > 1:
                raise ValueError('"z" argument must have only 1 non-degenerated dimension')
            else:
                dim = dim[0,0]
    elif z.ndim == 1:
        if dim is None:
            dim = np.argwhere(np.array(var_mask.shape) == z.size)
            if dim.size > 1:
                raise ValueError('Cannot identify vertical dimension')
            else:
                dim = dim[0,0]
                print('Note: dimension {:} considered as vertical, based on shape comparison'.format(dim))
    else:
        raise ValueError('"z" argument must be rank-1 or rank-3')

    # Check vertical compatibility and determine horizontal shape
    nz = z.size
    if var_mask.shape[dim] != nz:
        raise ValueError('vertical dimension of "var_mask" and "z" incompatible')

    horiz_shp = list(var_mask.shape)
    del(horiz_shp[dim])

    # Initialization
    bottom_idx = -1*np.ones(horiz_shp, dtype=int)
    idx = [slice(None), slice(None), slice(None)]

    # check z orientation
    if z[-1] < z[0]:
        zrange = range(nz-1, -1, -1)
    else:
        zrange = range(1, nz)

    # Loop on vertical dimension of "var_mask"
    for i in zrange:
        idx[dim] = i
        bottom_idx[~var_mask[tuple(idx)]] = i

    # Trick: add an extra element "0" and the end of z, so that if
    # bottom_idx==-1 (all elements are True) => z[bottom_idx] = 0
    z = z.reshape((nz,))
    z = np.concatenate((z, [0]))

    #<><><><><><><><><><>#
    depth = z[bottom_idx]
    #<><><><><><><><><><>#

    return depth


def get_var(*args, varname=None, raise_error=False):
    '''
    Try to get a variable "varname" (string) in one of the provided netCDF4 dataset
    (*args is all the datasets, provided as separated arguments).
    A None is got in the following cases:
      - dataset is None
      - varname is None
      - varname not found in dataset
    The first dataset that doesn't yield a None will be kept. If all do, return None.
    '''
    for dataset in args:
        if varname is None or dataset is None:
            var = None
        else:
            try:
                var = dataset[varname]
            except IndexError:
                var = None

        if var is not None:
            return var

    if raise_error:
        raise IndexError('Variable "'+str(varname)+'" cannot be found in any of the input datasets')



#################################################################
# MAIN FUNCTION: READ NETCDF FILES AND WRITE GEOCLIM INPUT FILE #
#################################################################


def write_oceanic_input(co2, input_files, latitude, z, temperature,
                        depth=None, time_dim=None, root=None, grid_file=None,
                        x_weight=None, y_weight=None, z_weight=None, horiz_weight=None, weight=None,
                        temp_outfile='GCM_oceanic_temp.dat',
                        surf_outfile='oce_surf.dat',
                        sedsurf_outfile='surf_sedi.dat',
                        vol_outfile='oce_vol.dat'):
    '''
    This function read the oceanic output and grid definition of a GCM,
    compute the volumes, surfaces and basin-average temperature of GEOCLIM
    basins, and save them in 4 GEOCLIM-readable ascii files, that can be
    used directly as GEOCLIM input: temperature file, box surface file,
    box sedimentary surface file, and box volume file.
    UP TO NOW, ONLY THE CREATION OF TEMPERATURE INPUT FILE IS IMPLEMENTED.

    mandatory input arguments:

        co2: list of CO2 levels (in increasing or decreasing order)
        input_files: list of names (string) of GCM output netCDF files,
                     corresponding to the CO2 levels (argument "co2")
        latitude: string. Name of the variable giving the latitude of
                  each point. That variable should be either in the
                  main input files, or in the (optional) grid file
                  (see argument "grid_file"). It can be 1D, or 2D
                  (if the horizontal grid is not longitude-latitude,
                  for instance)
        z: string. Name the variable giving the depth of each point
           of the 3D grid. That variable should be either in the main
           input files, or in the (optional) grid file (see argument
           "grid_file"). It must be 1D.
        temperature: string. Name the variable giving the temperature
                     of each point of the 3D grid. That variable should
                     be either in all the main input files (at each CO2
                     levels). Must be 3D or 4D if time-dependent. In that
                     case, it will be averaged along time-dimension.

    Optional input arguments:

        depth: string. Name the variable giving, at each "horizontal"
               grid point, the depth of ocean floor. That variable
               should be either in the main input files, or in the
               (optional) grid file (see argument "grid_file"). It
               must be 2D (horizontal).
               If not provided, the seafloor depth will be computed
               according to the valid points of temperature field.
        time_dim: string. Name of the time dimension if temperature
                  needs to be averaged on time. Note that the program
                  will by default try to identify it using common
                  time dimension names.
        root: string. Root directory where to the input files are (or
              more generally, common part to all file paths)
        grid_file: string. Path of an alternative file to get the grid
                   variables (latitude, z, depth and weightings)
        x_weight: string. Name of the variable giving the weighting of
                  the grid in the x direction (ie, width). Must be 1D
                  or 2D.
        y_weight: string. Name of the variable giving the weighting of
                  the grid in the y direction (ie, width). Must be 1D
                  or 2D.
        z_weight: string. Name of the variable giving the weighting of
                  the grid in the z direction (ie, thickness). Must be
                  1D.
        horiz_weight: string. Name of the variable giving the weighting of
                  the grid in the horizontal directions (x and y. ie, area).
                  Must be 2D.
        weight: string. Name of the variable giving the 3-D weighting of
                each point of the grid. Must be 3D.
        temp_outfile, surf_outfile, sedsurf_outfile, vol_outfile:
                name (ie, path) of the output files for (respectively)
                mean temperature, top surface, sedimentary surface, and
                volume, of the basins. See function definition for default
                names.
                Traditionally, surfaces and volume files are placed in
                INPUT/COMBINE/, whereas temperature file is placed in
                INPUT/YOUR_GCM_NAME/ (because it is CO2-dependent).

        Note, the priority order of the variables for the grid cell
        weighting (in case of conflict) is the following:
          weight > z_weight,horiz_weight > x_weight,y_weight
        In case of missing information for the weighting in one or several
        dimensions, A UNIFORM WEIGHTING WILL BE APPLIED.
    '''


    root = '' if root is None else root


    # CO2 is supposed to be in decreasing order:
    if co2[-1] > co2[0]:
        co2 = co2[::-1]
        input_files = input_files[::-1]


    # Grid information:
    # -----------------

    # Grid file (if not provided, get variables in the 1st main files)
    fgrid  = None if grid_file is None else nc.Dataset(grid_file)
    fgrid2 = nc.Dataset(root+input_files[0])

    xwg = get_var(fgrid, fgrid2, varname=x_weight)
    ywg = get_var(fgrid, fgrid2, varname=y_weight)
    zwg = get_var(fgrid, fgrid2, varname=z_weight)
    hwg = get_var(fgrid, fgrid2, varname=horiz_weight)
    wgh = get_var(fgrid, fgrid2, varname=weight)

    # Main variables (try the main file before grid file)
    lati = get_var(fgrid2, fgrid, varname=latitude, raise_error=True)
    zvar = get_var(fgrid2, fgrid, varname=z, raise_error=True)
    dpth = get_var(fgrid2, fgrid, varname=depth)

    # temperature of the 1st main file (to check existence, shape and dimensions)
    temp = fgrid2[temperature]

    # Check units (perform dummy conversion)
    _ = temperature_units.convert(0, temp.units)

    # Dimension checks:
    if temp.ndim == 3:
        idx_t = () # trick: empty tuple, so that temp.mean(idx_t) does nothing.
        temp_shape = temp.shape
        temp_dims = temp.dimensions
    elif temp.ndim == 4:
        if time_dim is None:
            idx_t = None
            for tname in ['time', 'time_counter', 'month']:
                if tname in temp.dimensions:
                    time_dim = tname
                    break

            if time_dim is None:
                raise ValueError('Cannot identify time dimension for temperature variable.')

        idx_t = temp.dimensions.index(time_dim)
        temp_shape = list(temp.shape)
        del(temp_shape[idx_t])
        temp_shape = tuple(temp_shape)
        temp_dims = list(temp.dimensions)
        del(temp_dims[idx_t])
        temp_dims = tuple(temp_dims)
                
    else:
        raise ValueError('temperature variable (argument "temperature") must be rank-3, or rank-4 if time-dependent.')

    if zvar.ndim != 1:
        raise ValueError('z variable (argument "z") must be rank-1')

    # Identify horizontal and vertical dimensions
    try:
        idx_z = temp_dims.index(zvar.dimensions[0])
    except ValueError:
        raise ValueError('Temperature not defined on the same dimension than z')

    if zvar.dimensions[0] in lati.dimensions:
        raise ValueError('latitude variable cannot be defined on z dimension')

    idx_xy = []
    for i in range(lati.ndim):
        try:
            idx_xy.append(temp_dims.index(lati.dimensions[i]))
        except ValueError:
            print('WARNING: temperature dimensions names do no match latitude ones')

    # create "nav_z" and "nav_lat"
    shp = [1,1,1]
    shp[idx_z] = zvar.shape[0]
    nav_z = length_units.convert(zvar[:], zvar.units).reshape(shp)

    # z must be positive
    if (nav_z<=0).all():
        nav_z = -nav_z

    shp = [1,1,1]
    for i0,i in enumerate(idx_xy):
        shp[i] = lati.shape[i0]

    nav_lat = latitude_units.convert(lati[:], lati.units).reshape(shp)

    # Get depth and create "nav_depth"
    if depth is None:
        print('Note: "depth" field automatically computed as the z-value of last valid point of temperature field on each vertical column')
        dpth = find_depth(temp[:].mean(idx_t).mask, nav_z)
    else:
        if dpth.ndim != 2:
            raise ValueError('seafloor depth variable (argument "depth") must be rank-2')
        else:
            for d in lati.dimensions:
                if d not in dpth.dimensions:
                    print('WARNING: latitude and seafloor depth variables dimensions does not have the same names')
            
            dpth = length_units.convert(dpth[:], dpth.units)

    shp = list(temp_shape)
    shp[idx_z] = 1
    nav_depth = dpth.reshape(shp)


    # Compute 3-D grid weighting
    # --------------------------

    if wgh is not None:
        if wgh.ndim != 3:
            raise ValueError('weight variable (argument "weight") must be rank-3')
        else:
            if wgh.shape != temp_shape:
                raise ValueError('temperature and weight variables (arguments "temperature" and "weight") must have the same shape')
            else:
                if temp_dims != wgh.dimensions: 
                    print('WARNING: temperature and weight variables dimensions does not have the same names')

                # + + + + + + #
                w = wgh[:,:,:]
                # + + + + + + #


    else:


        if hwg is not None:
            if hwg.ndim != 2:
                raise ValueError('horizontal weight variable (argument "horiz_weight") must be rank-2')
            else:
                for d in lati.dimensions:
                    if d not in hgh.dimensions:
                        print('WARNING: latitude and horizontal-weight variables dimensions does not have the same names')

                w_horiz = hwg[:,:]
                shp = list(temp_shape)
                shp[idx_z] = 1
                w_horiz = w_horiz.reshape(shp)

        else:

            if xwg is None and ywg is None:
                shp = list(temp_shape)
                shp[idx_z] = 1
                w_horiz = np.ones(shp, dtype=int)
                print('NOTE: uniform weighting assumed in horizontal dimensions (x and y)')
            else:

                idx_x = None
                idx_y = None

                if xwg is not None:
                    if xwg.ndim != 1:
                        raise ValueError('x weight variable (argument "x_weight") must be rank-1')
                    else:
                        try:
                            idx_x = temp_dims.index(xwg.dimensions[0])
                        except ValueError:
                            print('WARNING: x-weight variable dimension name does not correspond to temperature ones')

                        w_x = xwg[:]
                else:
                    w_x = np.ones(1, dtype=int)
                    print('NOTE: uniform weighting assumed in x dimension')

                if ywg is not None:
                    if ywg.ndim != 1:
                        raise ValueError('y weight variable (argument "u_weight") must be rank-1')
                    else:
                        try:
                            idx_y = temp_dims.index(ywg.dimensions[0])
                        except ValueError:
                            print('WARNING: y-weight variable dimension name does not correspond to temperature ones')

                        w_y = ywg[:]
                else:
                    w_y = np.ones(1, dtype=int)
                    print('NOTE: uniform weighting assumed in y dimension')

                # Check dimension order

                if idx_x is None:
                    idx_x = 0
                    while idx_x==idx_z or idx_x==idx_y:
                        idx_x += 1

                if idx_y is None:
                    idx_y = 0
                    while idx_y==idx_z or idx_y==idx_x:
                        idx_y += 1

                # Check eventual transposition
                if (w_x.size != 1 and w_x.size != temp_shape[idx_x] and w_x.size == temp_shape[idx_y]) or (w_y.size != 1 and w_y.size != temp_shape[idx_y] and w_y.size == temp_shape[idx_x]):
                    print('WARNING: transposition of x and y dimension assumed for computing horizontal weight')
                    idx_x,idx_y = idx_y,idx_x

                shp = [1,1,1]
                shp[idx_x] = temp_shape[idx_x]
                if w_x.size == 1:
                    w_x = np.ones(shp)
                else:
                    w_x = w_x.reshape(shp)

                shp = [1,1,1]
                shp[idx_y] = temp_shape[idx_y]
                if w_y.size == 1:
                    w_y = np.ones(shp)
                else:
                    w_y = w_y.reshape(shp)

                w_horiz = w_x * w_y


            if zwg is not None:
                if zwg.ndim != 1:
                    raise ValueError('z weight variable (argument "z_weight") must be rank-1')
                else:
                    if zwg.dimensions != zvar.dimensions:
                        print('WARNING: z-weight variable dimension name differs from z dimension name')

                    w_z = zwg[:]
            else:
                w_z = np.ones(zvar.shape, dtype=int)
                print('NOTE: uniform weighting assumed in z dimension')


            shp = [1,1,1]
            shp[idx_z] = w_z.size
            w_z = w_z.reshape(shp)

            # + + + + + + + #
            w = w_horiz * w_z
            # + + + + + + + #




    # Compute and export temperature average per basin -- loop on input files
    # -----------------------------------------------------------------------

    # Get mask of each basin of GEOCLIM
    Gmask = geoclim_basin_mask(nav_lat, nav_z, nav_depth)

    # Initialization
    nCO2 = len(co2)
    geoclim_input = np.zeros((nCO2,NBASIN+1), dtype=float)
    # Note: the first column of each row in GEOCLIM input file must be CO2 levels.
    #       the others columns are basin temperatures (in GEOCLIM order)

    ## Check: draw maps of masks
    ## from matplotlib import pyplot as plt
    ## temp = fgrid2[temperature]
    ## for i in range(NBASIN):
    ##     tempvar = temp[:].mean(idx_t)
    ##     print('Basin #{0:}, volume fraction: {1:f}'.format(i, w[np.logical_and(Gmask[i,:,:,:], ~tempvar.mask)].sum() / 
    ##                                                           w[~tempvar.mask].sum()))
    ##     plt.pcolormesh(np.logical_and(Gmask[i,:,:,:], ~(temp[:].mean(idx_t).mask)).any(idx_z))
    ##     plt.show()

    for k in range(nCO2):

        c = co2[k]
        fname = root+input_files[k]

        fin = nc.Dataset(fname)

        # Main variable
        temp = fin[temperature]

        # CO2 level
        geoclim_input[k,0] = c

        # basin-average temperature
        for i in range(NBASIN):
            tempvar = temp[:].mean(idx_t)
            geoclim_input[k,i+1] = (tempvar[Gmask[i,:,:,:]] * w[Gmask[i,:,:,:]]).sum() / w[np.logical_and(Gmask[i,:,:,:], ~tempvar.mask)].sum()

        geoclim_input[k,1:] = temperature_units.convert(geoclim_input[k,1:], temp.units)

        fin.close()
        #print(geoclim_input[k,:])


    # Write in file
    np.savetxt(temp_outfile, geoclim_input, fmt='%10.5f')


    if fgrid is not None:
        fgrid.close()



###################################################################################################



write_oceanic_input(co2=[560, 840, 1120, 1210, 1305, 1400, 1680, 1960, 2240, 4480, 8960],
                    input_files=['560ppm.nc', '840ppm.nc', '1120ppm.nc', '1210ppm.nc', '1305ppm.nc',
                                 '1400ppm.nc', '1680ppm.nc', '1960ppm.nc', '2240ppm.nc', '4480ppm.nc', '8960ppm.nc'],
                    root='/home/piermafrost/Downloads/150Ma/merg150_VegDef_AdjCSol_EccN_ocean_',
                    latitude='lat', z='lev', temperature='TEMP',
                    y_weight='lw', z_weight='thickness')

