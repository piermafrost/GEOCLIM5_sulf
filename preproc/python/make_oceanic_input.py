## Example of use:
##
## write_oceanic_temp(co2=[560, 840, 1120],
##                    input_files=['560ppm.nc', '840ppm.nc', '1120ppm.nc'],
##                    root='150Ma/merg150_VegDef_AdjCSol_EccN_ocean_',
##                    latitude='lat', z='lev', temperature='TEMP',
##                    y_weight='lw', z_weight='thickness')

import netCDF4 as nc
import numpy as np
from units import temperature_units, length_units, area_units, volume_units, latitude_units, velocity_units, flux_units, UnknownUnitError
from itertools import cycle
import os




###############################
# SPECIFIC GEOCLIM PARAMETERS #
###############################

# Oceanic basin definition
# ------------------------
NBASIN = 9 # not counting atmosphere box
H_EPICONT = 200  # (m) maximum ocean floor depth (bathymetry) of the epicontinental box
H_SURF    = 100  # (m) depth of surface boxes (ie, starting depth of thermocline or epicont deep box)
H_THERMO  = 1000 # (m) depth of thermocline boxes (ie, starting depth of deep box)
MIDLAT_RANGE = [-60, 60] # (degrees N): latitudinal range of inner ocean (ie, excluding high-latitude oceans)

# Modern total ocean volume and areas
# -----------------------------------
MODERN_VOLUME = 1.335e18 # m3
MODERN_AREA   = 3.619e14 # m2

# Arbitrary and and technical GEOCLIM parameters (must not be modified)
# ---------------------------------------------------------------------
ATM_VOL = 1.
ATM_AREA = 0.363e15
ATM_SEDI_AREA = 0.357e15
AREA_CONVERSION_FACTOR   = 1e-15 # areas expressed in 1e9 km2
VOLUME_CONVERSION_FACTOR = 1e-15 # volume expressed in 1e6 km3
FLUX_CONVERSION_FACTOR   = 1e-6  # water fluxes expressed in Sv (1e6 m3/s)


##########################
# MASK OF GEOCLIM BASINS #
##########################

def geoclim_basin_mask(nav_lat, nav_z, nav_bathy):
    '''
    Create a mask telling if points are in (T) or outside (F) each of the 9 GEOCLIM
    basins, for a given latitude array "nav_lat", local depth array "nav_z",
    and ocean floor depth (bathymetry) array "nav_bathy".
    the returned mask will rank-4, the 1st dimension being GEOCLIM basins, the last
    3 being the dimensions of the 3-D input fields (i.e., latitude, z and bathymetry)

    nav_lat, nav_z and nav_bathy MUST BE RANK-3, even if they are not defined on
    all dimensions (e.g., latitude and bathymetry should be defined on horizontal,
    dimensions, whereas z should be defined on the vertical dimension). Use
    degenerated (size-1) dimensions for those "extra" dimensions.
    The 3 dimensions must obviously correspond (ie, be in the same order) to the
    dimensions of 3D oceanic variables (e.g., temperature).

    The GEOCLIM basin definition must be the following:
      0. N high-lat, surface (incl. thermo)
      1. N high-lat, deep
      2. inner oce, surface
      3. inner oce, thermocline
      4. inner oce, deep
      5. epicont, surface
      6. epicont, deep
      7. S high-lat, surface (incl. thermo)
      8. S high-lat, deep
    '''

    i_surfbox = np.array([0, 2, 5, 7])

    shp = np.maximum(np.maximum(nav_lat.shape, nav_z.shape), nav_bathy.shape)

    mask = np.ones(np.concatenate(([NBASIN], shp)), dtype=bool)


    epicont_mask    = (nav_bathy <= H_EPICONT)
    surf_mask       = (nav_z <= H_SURF)
    surfthermo_mask = (nav_z <= H_THERMO)
    thermo_mask     = (nav_z > H_SURF)
    deep_mask       = (nav_z > H_THERMO)
    bottom_mask     = (nav_z <= nav_bathy)
    north_mask      = (nav_lat > MIDLAT_RANGE[1])
    south_mask      = (nav_lat < MIDLAT_RANGE[0])

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

    # Surface basins (excluding high latitude)
    for i in [2,5]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surf_mask)

    # Surface basins (high latitude)
    for i in [0,7]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surfthermo_mask)

    # Thermocline, and epicontinental deep basins
    for i in [3,6]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], thermo_mask)

    # Thermocline basin
    for i in [3]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surfthermo_mask)

    # Other ocean deep basin
    for i in [1,4,8]:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], deep_mask)

    ## + epicontinental (remove points below seafloor)
    #for i in [1,4,5,6,8]:
    #    mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], bottom_mask)

    return mask, i_surfbox


##########################################
# FUNCTIONS SPECIFIC TO PARTICULAR GRIDS #
##########################################

def remove_orca_north_fold(x):
    '''
    Put 0 in points outside the irregular inner domain mask
    at North fold boundary of ORCA grid
    '''
    ihalf = x.shape[-1]//2
    if x.ndim == 2:
        x[-1, :] = 0
        x[-2, :ihalf+1] = 0
    elif x.ndim == 3:
        x[:, -1, :] = 0
        x[:, -2, :ihalf+1] = 0


#######################
# AUXILIARY FUNCTIONS #
#######################

def find_bathy(var_mask: np.ndarray, z: np.ndarray, invert_zaxis=False):
    '''
    Return a 2D (horizontal) array giving the seafloor depth (bathymetry) of each point,
    and the 2D array of vertical indices at which the seafloor is reached.
    The bathymetry is determined as follows:
      1. Find the index the last 'False' of 3D boolean array "var_mask" in the vertical
         (z) dimension (yield a 2D "horizontal" index array `idx`).
      2. Pick the value of "z" for each "bottom index", ie, `z[i,j,idx[i,j]]` for all i,j
         If all the value of var_mask in one column are 'True', the depth of that point
         will be '0'
    "var_mask" must be a rank-3 boolean array. "z" must be a rank-3 array, and can have
    degenerated (size-1) dimensions. The vertical dimension must be the 1st one (C-indexing),
    must be positive downward, and must be in increasing order (i.e., from surface to
    bottom), unless specified "invert_zaxis=True"
    '''

    # order of z axis
    order = slice(None, None, -1) if invert_zaxis else slice(None)

    if var_mask.ndim != 3 or z.ndim != 3:
        raise ValueError('"var_mask" and "z" must be rank-3')

    # Check vertical compatibility and determine horizontal shape
    if not (z.shape[0]==var_mask.shape[0] and \
            np.logical_or(np.array(z.shape)==1, np.array(var_mask.shape)==np.array(z.shape)).all()):
        raise ValueError('vertical dimension of "var_mask" and "z" incompatible')

    nz = z.shape[0]
    horiz_shp = var_mask.shape[1:]

    # Initialization
    if invert_zaxis:
        bottom_idx = -1*np.ones(horiz_shp, dtype='int16')
    else:
        bottom_idx = nz*np.ones(horiz_shp, dtype='int16')

    # Loop on vertical dimension of "var_mask"
    for i in range(nz)[order]:
        bottom_idx[~var_mask[i,:,:]] = i

    # Trick: add an extra element "0" and the end of z, so that if
    # bottom_idx=="last index" (all elements are True) => z[bottom_idx] = 0
    if invert_zaxis:
        z = np.concatenate((np.zeros((1,)+z.shape[1:], z.dtype), z), axis=0)
        bottom_idx += 1
    else:
        z = np.concatenate((z, np.zeros((1,)+z.shape[1:], z.dtype)), axis=0)

    #<><><><><><><><><><><><><><><>#
    ii, jj = np.indices(z.shape[1:])
    bathy = z[bottom_idx, ii, jj]
    #<><><><><><><><><><><><><><><>#

    if invert_zaxis:
        bottom_idx -= 1
        # Replace "last index" values by "0"
        bottom_idx[bottom_idx==-1] = 0
    else:
        # Replace "last index" values by "-1"
        bottom_idx[bottom_idx==nz] = -1

    return bathy, bottom_idx

# ========== #

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

# ========== #

def get_axis_weighting(axis: str, wghvar, ref_dims, ref_shp, axis_bnds_var=None, axis_var=None, t_dim=None):
    '''
    Get the weighting of one axis (x, y or z), check that the dimensions and size
    match the 3D reference shape.
    input argument:
      - "axis": name of the axis (string)
      - "wghvar": netCDF4 variable (in netCDF4 Dataset) of the weighting.
        > if wghvar is None: create uniform weighting 
      - "ref_dims": tuple of dimensions of the reference 3D variable
      - "ref_shp": shape (tuple) of the reference 3D variable arrays
    optional input arguments:
      - "axis_bnds_var": netCDF4 variable (in netCDF4 Dataset) giving the axis 
          bounds (must be 2D, with 2nd dimension, size-2, being "bounds").
          Get axis weighting as the difference (np.diff) of axis bounds in case
          "wghvar" is None.
      - "axis_var": netCDF4 variable (in netCDF4 Dataset) giving the axis
          coordinates. Used to get "axis_bnds_var" units in case "axis_bnds_var"
          doesn't have a "units" attribute. 
      - "t_dim": potential time dimension in case the "wghvar" is 4D, and will
          be averaged on the time dimension.
    returns:
      - 3D array (with shape "ref_shp") of the axis weighting
      - boolean, indicated whether axis weighting has units (expect "meter")
    '''

    if wghvar is not None:

        shp = list(ref_shp)

        if wghvar.ndim == 3:
            wgh = wghvar[:,:,:]

        elif wghvar.ndim == 2:
            if wghvar.dimensions[0]==ref_dims[1] and wghvar.dimensions[1]==ref_dims[2]: # weighting dimensions are {x,y}
                shp[0] = 1 # remove z dimension
            elif wghvar.dimensions[0]==ref_dims[0] and wghvar.dimensions[1]==ref_dims[2]: # weighting dimensions are {x,z}
                shp[1] = 1 # remove y dimension
            elif wghvar.dimensions[0]==ref_dims[0] and wghvar.dimensions[1]==ref_dims[1]: # weighting dimensions are {y,z}
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_weight" variable (do not match "temperature" variable)')

            wgh = wghvar[:,:].reshape(shp)

        elif wghvar.ndim == 1:
            if wghvar.dimensions[0]==ref_dims[2]: # weighting dimension is {x}
                shp[0] = 1 # remove z dimension
                shp[1] = 1 # remove y dimension
            elif wghvar.dimensions[0]==ref_dims[1]: # weighting dimension is {y}
                shp[0] = 1 # remove z dimension
                shp[2] = 1 # remove x dimension
            elif wghvar.dimensions[0]==ref_dims[0]: # weighting dimension is {z}
                shp[1] = 1 # remove y dimension
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_weight" variable (do not match "temperature" variable)')

            wgh = wghvar[:].reshape(shp)

        else:
            if wghvar.ndim == 4 and t_dim is not None:
                wgh = wghvar[:,:,:,:].mean(t_dim)
            else:
                raise ValueError('Illegal number of dimensions for variable "'+axis+'_weight"')

        try:
            wgh = length_units.convert(wgh, wghvar.units)
            has_units = True
        except (AttributeError, UnknownUnitError):
            has_units = False

    elif axis_bnds_var is not None:
        if axis_bnds_var.ndim != 2:
            raise ValueError('Cannot handle "'+axis+'_bounds" variable that is not 2D')
        elif axis_bnds_var.shape[1] != 2:
            raise ValueError('2nd dimension of '+axis+' bounds variable (argument "'+axis+'_bounds") must be size-2')
        else:
            wgh = np.abs(np.diff(axis_bnds_var[:], axis=1))
            try:
                wgh = length_units.convert(wgh, axis_bnds_var.units)
                has_units = True
            except AttributeError:
                if axis_var is None:
                    has_units = False
                else:
                    try:
                        wgh = length_units.convert(wgh, axis_var.units)
                        has_units = True
                        print('WARNING: units of '+axis+' bounds variable not found. ASSUME SAME UNIT AS '+axis+' VARIABLE.')
                    except UnknownUnitError:
                        has_units = False

            except UnknownUnitError:
                has_units = False

            shp = list(ref_shp)
            if axis_bnds_var.dimensions[0]==ref_dims[2]: # weighting dimension is {x}
                shp[0] = 1 # remove z dimension
                shp[1] = 1 # remove y dimension
            elif axis_bnds_var.dimensions[0]==ref_dims[1]: # weighting dimension is {y}
                shp[0] = 1 # remove z dimension
                shp[2] = 1 # remove x dimension
            elif axis_bnds_var.dimensions[0]==ref_dims[0]: # weighting dimension is {z}
                shp[1] = 1 # remove y dimension
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_bounds" variable (do not match "temperature" variable)')

            wgh = wgh.reshape(shp)

    else:
        wgh = 1
        print('WARNING: uniform weighting assumed in '+axis+' dimension')
        has_units = False

    return wgh*np.ones(ref_shp), has_units

# ========== #

def average(*args):
    '''
    Return the unweighted average of all input arguments (that must be
    numpy.ma.core.MaskedArray of shame shape), removing masked value from
    the average. Where all arrays are masked, the average is 0.
    '''
    ave = np.zeros(args[0].shape, dtype=float)
    cnt = np.zeros(args[0].shape, dtype='int8')
    msk = np.zeros(args[0].shape, dtype=bool)
    for x in args:
        #x = np.ma.array(x)
        ave[~x.mask] += x.data[~x.mask]
        cnt[~x.mask] += 1
        msk[~x.mask] = True

    ave[msk] /= cnt[msk]
    return ave

# ========== #

def minimum(*args):
    '''
    Return the point-by-point minimum of all input arguments (that must be
    numpy.ma.core.MaskedArray of shame shape), removing masked value from
    the list to compute the miminum.
    Where all arrays are masked, the minimum is 0.
    '''
    mini = np.zeros(args[0].shape, dtype=float)
    msk  = np.zeros(args[0].shape, dtype=bool)
    for x in args:
        #x = np.ma.array(x)
        msk1 = np.logical_and(~x.mask, ~msk)
        msk2 = np.logical_and(~x.mask, msk)
        mini[msk1] = x.data[msk1]
        mini[msk2] = np.minimum(mini[msk2], x.data[msk2])
        msk[msk1] = True

    return mini



#################################################################
# MAIN FUNCTION: READ NETCDF FILES AND WRITE GEOCLIM INPUT FILE #
#################################################################


def write_oceanic_input(co2, input_files, latitude, z, temperature,
                        bathy=None, time_dim=None, root='', grid_file=None,
                        x_weight=None, y_weight=None, z_weight=None, z_bounds=None, horiz_weight=None, weight=None,
                        input_u_files=None, input_v_files=None, input_w_files=None,              # !!
                        u=None, v=None, w=None, u_flx=None, v_flx=None, w_flx=None,              # !!
                        inverted_w=False, periodic_x=True, x_overlap=0, special_wrap=None,       # !! NEW ARGUMENTS
                        outdir='./',                                                             # !!
                        temp_outfile='GCM_oceanic_temp.dat',
                        surf_outfile='oce_surf.dat',
                        sedsurf_outfile='surf_sedi.dat',
                        vol_outfile='oce_vol.dat',
                        watflux_outfile='exchange_2.dat'):
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

        bathy: string. Name the variable giving, at each "horizontal"
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
                   variables (latitude, z, bathy and weightings)
        x_weight: string. Name of the variable giving the weighting of
                  the grid in the x direction (ie, width). Must be 1D
                  or 2D.
        y_weight: string. Name of the variable giving the weighting of
                  the grid in the y direction (ie, width). Must be 1D
                  or 2D.
        z_weight: string. Name of the variable giving the weighting of
                  the grid in the z direction (ie, thickness). Must be
                  1D.
        z_bounds: string. Alternative to "z_weight", name of the variable
                  giving the upper and lower bounds of the vertical levels.
                  z-weighting will then be computed as diff(z_bounds).
                  Must be 2D, the 2nd dimension being the bounds (size-2).
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

        Notes:

        the priority order of the variables for the grid cell weighting
        (in case of conflict) is the following:
          weight > z_weight,horiz_weight > x_weight,y_weight
        In case of missing information for the weighting in one or several
        dimensions, A UNIFORM WEIGHTING WILL BE APPLIED.

        Output area unit is Gkm2 (1e15 m2)
        Output volume unit is Mkm3 (1e15 m3)

        The volume of the atmospheric box is expected to be 1e-15 (so that
        is is "1" after the conversion in m3)

        Example:

        write_oceanic_temp(co2=[560, 840, 1120],
                           input_files=['560ppm.nc', '840ppm.nc', '1120ppm.nc'],
                           root='150Ma/merg150_VegDef_AdjCSol_EccN_ocean_',
                           latitude='lat', z='lev', temperature='TEMP',
                           y_weight='lw', z_weight='thickness')
    '''


    # Output directory
    if outdir[-1]!='/':
        outdir = outdir+'/'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        print('Created directory "'+outdir+'"')


    # CO2 array (for GEOCLIM without climatic parameters)
    # ---------------------------------------------------

    if co2 is None:

        nclim = len(input_files)
        co2 = cycle([-1e36])

    else:

        nclim = len(co2)

        # CO2 is supposed to be in decreasing order:
        if co2[-1] > co2[0]:
            co2 = co2[::-1]
            input_files   = input_files[::-1]
            input_u_files = input_u_files[::-1]
            input_v_files = input_v_files[::-1]
            input_w_files = input_w_files[::-1]


    # Number of input files
    # ---------------------

    if len(input_files) != nclim:
        raise ValueError('Number of input files inconsistent with CO2 axis')

    if input_u_files is None:
        input_u_files = cycle([None])
    elif len(input_u_files) != nclim:
        raise ValueError('Number of input "u" files inconsistent with number of input "Temp" files')

    if input_v_files is None:
        input_v_files = cycle([None])
    elif len(input_v_files) != nclim:
        raise ValueError('Number of input "v" files inconsistent with number of input "Temp" files')

    if input_w_files is None:
        input_w_files = cycle([None])
    elif len(input_w_files) != nclim:
        raise ValueError('Number of input "w" files inconsistent with number of input "Temp" files')

    all_inputs = (input_files, input_u_files, input_v_files, input_w_files)


    # Grid information:
    # -----------------

    # Grid file (if not provided, get variables in the 1st main files)
    fgrid  = None if grid_file is None else nc.Dataset(grid_file)
    fgrid2 = nc.Dataset(root+input_files[0])

    xwg = get_var(fgrid, fgrid2, varname=x_weight)
    ywg = get_var(fgrid, fgrid2, varname=y_weight)
    zwg = get_var(fgrid, fgrid2, varname=z_weight)
    zbd = get_var(fgrid, fgrid2, varname=z_bounds)
    hwg = get_var(fgrid, fgrid2, varname=horiz_weight)
    wgh = get_var(fgrid, fgrid2, varname=weight)

    # Main variables (try the main file before grid file)
    lati = get_var(fgrid2, fgrid, varname=latitude, raise_error=True)
    zvar = get_var(fgrid2, fgrid, varname=z, raise_error=True)
    baty = get_var(fgrid2, fgrid, varname=bathy)

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


    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    #@ Expected dimensions (time excluded): @#
    #@  0: z                                @#
    #@  1: y                                @#
    #@  2: x                                @#
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


    # Special grid boundary conditions (x-overlapping, ...)
    # -----------------------------------------------------

    if special_wrap == 'ORCA':
        ix0 = 1
        ix1 = temp_shape[2]-1
        x_overlap = 2

    else:
        ix0 = 0
        ix1 = temp_shape[2] - x_overlap
        if special_wrap is not None:
            print('WARNING: unknown special wrapping case "'+str(special_wrap)+'".')
            print('         Will be ignored.')


    horiz_shp = temp_shape[1:]


    # create "nav_z"
    # --------------
    shp = list(temp_shape)
    if zvar.ndim == 3:
        nav_z = zvar[:,:,:]

    elif zvar.ndim == 2:
        if zvar.dimensions[1] == temp_dims[2]: # "z" dimensions are {x,z}
            shp[1] = 1 # remove y dimension
        elif zvar.dimensions[1] == temp_dims[1]: # "z" dimensions are {y,z}
            shp[2] = 1 # remove x dimension
        else:
            raise ValueError('Cannot identify dimensions of "z" variable (do not match "temperature" variable)')

        nav_z = zvar[:,:].reshape(shp)

    elif zvar.ndim == 1:
        shp[1] = 1 # remove y dimension
        shp[2] = 1 # remove x dimension
        nav_z = zvar[:].reshape(shp)

    else:
        raise ValueError('Illegal number of dimensions for variable "z"')

    # z must be positive
    if (nav_z<=0).all():
        nav_z = -nav_z

    # check z ordering:
    reverse_z = (nav_z[0,:,:] > nav_z[1,:,:]).any()
    isurf = -1 if reverse_z else 0

    # Units conversion
    nav_z = length_units.convert(nav_z, zvar.units)

    # create "nav_lat"
    # ----------------
    shp = list(temp_shape)
    if lati.ndim == 3:
        raise ValueError('Cannot handle 3-dimensionnal latitude field. Must be 1D or 2D')
    elif lati.ndim == 2:
        shp[0] = 1 # no z dimension
    elif lati.ndim == 1:
        shp[0] = 1 # no z
        shp[2] = 1 # nor x dimension
    else:
        raise ValueError('Illegal number of dimensions for variable "latitude"')

    nav_lat = latitude_units.convert(lati[:], lati.units).reshape(shp)

    # Get depth and create "nav_bathy"
    # --------------------------------
    if bathy is None:
        baty, bottom_idx = find_bathy(temp[:].mean(idx_t).mask, nav_z, invert_zaxis=reverse_z)
        print('NOTE: "bathy" field (bathymetry) automatically computed as the z-value of')
        print('      the last valid point of temperature field on each vertical column')
    else:
        _, bottom_idx = find_bathy(temp[:].mean(idx_t).mask, nav_z, invert_zaxis=reverse_z)
        if baty.ndim != 2:
            raise ValueError('seafloor depth variable (argument "bathy") must be rank-2')
        else:
            for d in lati.dimensions:
                if d not in baty.dimensions:
                    print('WARNING: latitude and seafloor depth variables dimensions does not have the same names')
            
            baty = length_units.convert(baty[:], baty.units)

    shp = list(temp_shape)
    shp[0] = 1 # no z dimension
    nav_bathy = baty.reshape(shp)


    # Volumetric (3-D) and/or area (horizontal) grid weighting
    # --------------------------------------------------------

    # Volumetric
    # ----------
    if wgh is not None:
        if wgh.ndim != 3:
            raise ValueError('weight variable (argument "weight") must be rank-3')
        else:
            if wgh.shape != temp_shape:
                raise ValueError('temperature and weight variables (arguments "temperature" and "weight") must have the same shape')
            else:
                if temp_dims != wgh.dimensions: 
                    print('WARNING: temperature and weight variables dimensions does not have the same names')

                w_vol = wgh[:,:,:]
                try:
                    w_vol = volume_units.convert(w_vol, wgh.units)
                    volunits = True
                except (AttributeError, UnknownUnitError):
                    volunits = False

    else:
        w_vol = None

    # Area
    # ----
    if hwg is None:
        w_horiz = None
    else:
        if hwg.ndim != 2:
            raise ValueError('horizontal weight variable (argument "horiz_weight") must be rank-2')
        else:
            for d in lati.dimensions:
                if d not in hwg.dimensions:
                    print('WARNING: latitude and horizontal-weight variables dimensions does not have the same names')

            w_horiz = hwg[:,:]
            try:
                w_horiz = area_units.convert(w_horiz, hwg.units)
                xyunits = True
            except (AttributeError, UnknownUnitError):
                xyunits = False


    # x, y and z axis weighting
    # -------------------------

    if xwg is None and ywg is None:
        if w_horiz is None:
            w_horiz = np.ones(horiz_shp)
            xyunits = False
            if w_vol is None:
                print('WARNING: uniform weighting assumed in horizontal dimensions (x and y)')

    # x weighting
    w_x, xunits = get_axis_weighting('x', xwg, temp_dims, temp_shape, t_dim=idx_t)

    # y weighting
    w_y, yunits = get_axis_weighting('y', ywg, temp_dims, temp_shape, t_dim=idx_t)

    # x weighting
    w_z, zunits = get_axis_weighting('z', zwg, temp_dims, temp_shape, axis_bnds_var=zbd, axis_var=zvar, t_dim=idx_t)


    # Compute areal (horizontal) weighting if it wasn't loaded
    if w_horiz is None:
        if xwg is None and ywg is None:
            xyunits = False
            if w_vol is None:
                print('WARNING: uniform weighting assumed in horizontal dimensions (x and y)')
            else:
                # consider first "surface" volumetric weighting
                w_horiz = w_vol[isurf,:,:]
                
        else:
            w_horiz = w_x[isurf,:,:] * w_y[isurf,:,:]
            xyunits = (xunits and yunits)

        if not xyunits:
            print('WARNING: x and y units not understood. Areas will be computed as a fraction of modern oceanic areas')


    # Compute volumetric weighting if it wasn't loaded
    if w_vol is None:
        w_vol = w_horiz.reshape((1,)+horiz_shp) * w_z
        volunits = (xyunits and zunits)
        if not volunits:
            print('WARNING: horizontal or vertical units not understood. Volume will be computed as a fraction of modern oceanic areas')


    # Special domain definition
    # + + + + + + + + + + + + +
    if special_wrap == 'ORCA':
        # => irregular inner domain mask at North fold boundary in ORCA grid
        remove_orca_north_fold(w_horiz)
        remove_orca_north_fold(w_vol)
        remove_orca_north_fold(w_x)
        remove_orca_north_fold(w_y)
        remove_orca_north_fold(w_z)


    # Compute and export temperature average per basin -- loop on input files
    # -----------------------------------------------------------------------

    # Get mask of each basin of GEOCLIM
    Gmask, i_surfbox = geoclim_basin_mask(nav_lat, nav_z, nav_bathy)

    # Initialization
    geoclim_input = np.zeros((nclim,NBASIN+1), dtype=float)
    # Note: the first column of each row in GEOCLIM input file must be CO2 levels.
    #       the others columns are basin temperatures (in GEOCLIM order)

    ## Check: draw maps of masks
    ## <+++++++++++++++++++++++++++++++++++++++> #
    #from matplotlib import pyplot as plt
    #for k in range(NBASIN):
    #    print('Basin #{0:}, volume fraction: {1:f}'.format(k, w_vol[Gmask[k,:,:,:]].sum()/w_vol.sum()))
    #    fig, ax = plt.subplots(nrows=3, figsize=(6,8))
    #    pid = ax[0].pcolormesh(np.count_nonzero(Gmask[k,:,:,:], axis=0))
    #    ax[0].set_xlabel('x')
    #    ax[0].set_ylabel('y')
    #    plt.colorbar(pid)
    #    pid = ax[1].pcolormesh(np.count_nonzero(Gmask[k,:,:,:], axis=2))
    #    ax[1].set_xlabel('y')
    #    ax[1].set_ylabel('z')
    #    plt.colorbar(pid)
    #    pid = ax[2].pcolormesh(np.count_nonzero(Gmask[k,:,:,:], axis=1))
    #    ax[2].set_xlabel('x')
    #    ax[2].set_ylabel('z')
    #    plt.colorbar(pid)
    #    fig.suptitle('Basin #{0:}'.format(k))
    #    if not reverse_z:
    #        ax[1].invert_yaxis()
    #        ax[2].invert_yaxis()
    #plt.show()
    ## <+++++++++++++++++++++++++++++++++++++++> #

    for k,(inputfile,c) in enumerate(zip(input_files, co2)):

        fname = root+inputfile

        fin = nc.Dataset(fname)

        # Main variable
        temp = fin[temperature]

        # CO2 level
        geoclim_input[k,0] = c

        # basin-average temperature
        for i in range(NBASIN):
            tempvar = temp[:].mean(idx_t)
            geoclim_input[k,i+1] = (tempvar[:,:,ix0:ix1][Gmask[i,:,:,ix0:ix1]] * w_vol[:,:,ix0:ix1][Gmask[i,:,:,ix0:ix1]]).sum() / \
                                   w_vol[:,:,ix0:ix1][np.logical_and(Gmask[i,:,:,ix0:ix1], ~tempvar.mask[:,:,ix0:ix1])].sum()

        geoclim_input[k,1:] = temperature_units.convert(geoclim_input[k,1:], temp.units)

        fin.close()


    # Write in file
    np.savetxt(outdir+temp_outfile, geoclim_input, fmt='%10.5f')


    # Compute and export water fluxes
    # -------------------------------


    if ((u is None and u_flx is None) or (v is None and v_flx is None) or (w is None and w_flx is None)):
        print('Missing "u", "u_flx", "v", "v_flx", "w" or "w_flx" fields. Exchange fluxes not computed.')

    elif ((u_flx is None and not (yunits and xunits)) or \
          (v_flx is None and not (xunits and zunits)) or \
          (w_flx is None and not (xunits and yunits))):
        print('Warning: cannot compute water fluxes: missing units for x, y or z axis weigthing')

    else:
        fout = open(outdir+watflux_outfile, mode='w')

        # Finite volume approach:
        # compute areas of the surfaces between cells (in the perpendiular direction)
        #    -> "y-z" area between cell x=i and cell x=i+1 is ("y-z" area[cell i]  +  "y-z" area[cell i+1])/2
        # except that for vertical dimension, the minimum "thickness" between cell i and cell i+1 is considered:
        #    vertical cell thickness is expected to be uniform along x and y dimension, except for bottom cells, that intercept
        #    sea-floor with different depths. -> take the minimum thickness between 2 adjacent cell as the exchange surface

        if u_flx is None: # average/minimum in the x direction
            w_yz = average(w_y[:,:,ix0:ix1-1], w_y[:,:,ix0+1:ix1]) * minimum(w_z[:,:,ix0:ix1-1], w_z[:,:,ix0+1:ix1])
            if periodic_x:
                w_yz_edge = average(w_y[:,:,ix1-1], w_y[:,:,ix0]) * minimum(w_z[:,:,ix1-1], w_z[:,:,ix0])

        if v_flx is None: # average/minimum in the y direction
            w_xz = average(w_x[:,:-1,ix0:ix1], w_x[:,1:,ix0:ix1]) * minimum(w_z[:,:-1,ix0:ix1], w_z[:,1:,ix0:ix1])

        if w_flx is None: # average in the z direction
            w_xy = average(w_x[:-1,:,ix0:ix1], w_x[1:,:,ix0:ix1]) * average(w_y[:-1,:,ix0:ix1], w_y[1:,:,ix0:ix1])

        units = lambda flx: velocity_units if flx is None else flux_units


        from matplotlib import pyplot as plt
        # Loop on all u,v,w input files => compute exchange matrix for all climate states
        # -------------------------------------------------------------------------------
        for allfiles in zip(*all_inputs):

            allf = [None if fname is None else nc.Dataset(root+fname) for fname in allfiles]

            u_var = get_var(*allf, varname=(u if u_flx is None else u_flx), raise_error=True)
            v_var = get_var(*allf, varname=(v if v_flx is None else v_flx), raise_error=True)
            w_var = get_var(*allf, varname=(w if w_flx is None else w_flx), raise_error=True)

            u_field = units(u_flx).convert(u_var[:].mean(idx_t), u_var.units)
            v_field = units(v_flx).convert(v_var[:].mean(idx_t), v_var.units)
            w_field = units(w_flx).convert(w_var[:].mean(idx_t), w_var.units)

            if special_wrap == 'ORCA':
                if u_flx is not None:
                    remove_orca_north_fold(u_field)

                if v_flx is not None:
                    remove_orca_north_fold(v_field)

                if w_flx is not None:
                    remove_orca_north_fold(w_field)


            # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
            # Matrix of fluxes: X_mat[i,j] = flux FROM box i TO box j #
            # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
            X_mat = np.zeros((NBASIN,NBASIN), dtype=float)

            for i in range(NBASIN):
                for j in list(range(i))+list(range(i+1,NBASIN)):

                    uij, vij, wij = (), (), ()

                    # u-fluxes
                    if periodic_x:
                        msk_uij = np.logical_and(Gmask[i,:,:,ix1-1], Gmask[j,:,:,ix0])
                        msk_uij = np.logical_and(msk_uij, ~u_field[:,:,ix1-1].mask)
                        if msk_uij.any():
                            uij = u_field[:,:,ix1-1][msk_uij].data
                            if u_flx is None:
                                uij *= w_yz_edge[msk_uij]

                            pos = (uij >= 0)
                            X_mat[i,j] += uij[pos].sum()
                            X_mat[j,i] -= uij[~pos].sum()

                    msk_uij = np.logical_and(Gmask[i,:,:,ix0:ix1-1], Gmask[j,:,:,ix0+1:ix1])
                    msk_uij = np.logical_and(msk_uij, ~u_field[:,:,ix0:ix1-1].mask)
                    if msk_uij.any():
                        uij = u_field[:,:,ix0:ix1-1][msk_uij].data
                        if u_flx is None:
                            uij *= w_yz[msk_uij]

                        pos = (uij >= 0)
                        X_mat[i,j] += uij[pos].sum()
                        X_mat[j,i] -= uij[~pos].sum()

                    # v-fluxes
                    msk_vij = np.logical_and(Gmask[i,:,:-1,ix0:ix1], Gmask[j,:,1:,ix0:ix1])
                    msk_vij = np.logical_and(msk_vij, ~v_field[:,:-1,ix0:ix1].mask)
                    if msk_vij.any():
                        vij = v_field[:,:-1,ix0:ix1][msk_vij].data
                        if v_flx is None:
                            vij *= w_xz[msk_vij]

                        pos = (vij > 0)
                        X_mat[i,j] += vij[pos].sum()
                        X_mat[j,i] -= vij[~pos].sum()

                    # w-fluxes
                    if inverted_w:
                        msk_wij = np.logical_and(Gmask[i,1:,:,ix0:ix1], Gmask[j,:-1,:,ix0:ix1])
                        msk_wij = np.logical_and(msk_wij, ~w_field[1:,:,ix0:ix1].mask)
                    else:
                        msk_wij = np.logical_and(Gmask[i,:-1,:,:], Gmask[j,1:,:,ix0:ix1])
                        msk_wij = np.logical_and(msk_wij, ~w_field[:-1,:,ix0:ix1].mask)

                    if msk_wij.any():

                        if inverted_w:
                            wij = w_field[1:,:,ix0:ix1][msk_wij].data
                        else:
                            wij = w_field[:-1,:,ix0:ix1][msk_wij].data

                        if w_flx is None:
                            wij *= w_xy[msk_wij]

                        pos = (wij >= 0)
                        X_mat[i,j] += wij[pos].sum()
                        X_mat[j,i] -= wij[~pos].sum()

                    ## <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #
                    #if msk_uij.any() or msk_vij.any() or msk_wij.any():
                    #    fig, ax = plt.subplots(3, 3, figsize=(15,8))
                    #    for k,(uvw,msk,txt) in enumerate(zip([uij, vij, wij], [msk_uij, msk_vij, msk_wij], 'UVW')):
                    #        x = np.zeros(msk.shape)
                    #        x[msk] = np.abs(uvw)
                    #        #pid = ax[0,k].pcolormesh(np.count_nonzero(msk, axis=0))
                    #        pid = ax[0,k].pcolormesh(np.ma.masked_values(x.max(axis=0), 0))
                    #        ax[0,k].set_xlabel('x')
                    #        ax[0,k].set_ylabel('y')
                    #        plt.colorbar(pid)
                    #        #pid = ax[1,k].pcolormesh(np.count_nonzero(msk, axis=2))
                    #        pid = ax[1,k].pcolormesh(np.ma.masked_values(x.max(axis=2), 0))
                    #        ax[1,k].set_xlabel('y')
                    #        ax[1,k].set_ylabel('z')
                    #        plt.colorbar(pid)
                    #        #pid = ax[2,k].pcolormesh(np.count_nonzero(msk, axis=1))
                    #        pid = ax[2,k].pcolormesh(np.ma.masked_values(x.max(axis=1), 0))
                    #        ax[2,k].set_xlabel('x')
                    #        ax[2,k].set_ylabel('z')
                    #        plt.colorbar(pid)
                    #        if not reverse_z:
                    #            ax[1,k].invert_yaxis()
                    #            ax[2,k].invert_yaxis()
                    #        ax[0,k].set_title(txt)
                    #    fig.suptitle('Basin {0:} => Basin {1:}'.format(i,j))
                    #    plt.show()
                    ## <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #

            np.savetxt(fout, FLUX_CONVERSION_FACTOR*X_mat, fmt='%.2f', delimiter='\t') # '%.7e'
            fout.write('\n')

        fout.close()


    if fgrid is not None:
        fgrid.close()


    # >                                                                                < #
    # > Following input files depend only on oceanic basins definition, not on climate < #
    # > ============================================================================== < #


    # Update basin mask to remove missing points of last loaded temperature
    for k in range(NBASIN):
        Gmask[k,:,:,:][tempvar.mask] = False


    # Compute and export box areas
    # ----------------------------

    area_top = np.zeros((NBASIN+1,), dtype=w_horiz.dtype)
    area_sed = np.zeros((NBASIN+1,), dtype=w_horiz.dtype)

    # "area_top" is the horizontal area at the top of boxes
    for k in range(NBASIN):
        area_top[k] = w_horiz[:,ix0:ix1][Gmask[k,:,:,ix0:ix1].any(0)].sum()

    # "area_sed" is the area where boxes intercep seafloor
    jj, ii = np.indices(np.array(horiz_shp) - np.array([0,x_overlap]))
    for k in range(NBASIN):
        area_sed[k] = w_horiz[:,ix0:ix1][Gmask[k,:,:,ix0:ix1][bottom_idx[:,ix0:ix1], jj, ii]].sum()

    if not xyunits:
        # normalize areas and scale them to modern total ocean volume
        print('areas scaled to modern ocean area')
        area_top *= MODERN_AREA/area_top[i_surfbox].sum()
        area_sed *= MODERN_AREA/area_sed.sum()

    # atmosphere box
    area_top[-1] = ATM_AREA
    area_sed[-1] = ATM_SEDI_AREA

    # Conversion to GEOCLIM units
    area_top *= AREA_CONVERSION_FACTOR
    area_sed *= AREA_CONVERSION_FACTOR

    # Write in file
    np.savetxt(outdir+surf_outfile,    area_top, fmt='%.12e', delimiter='\n')
    np.savetxt(outdir+sedsurf_outfile, area_sed, fmt='%.12e', delimiter='\n')


    # Compute and export box volumes
    # ------------------------------

    box_vol = np.zeros((NBASIN+1), dtype=w_vol.dtype)
    for k in range(NBASIN):
        box_vol[k] = w_vol[:,:,ix0:ix1][Gmask[k,:,:,ix0:ix1]].sum()

    if not volunits:
        # normalize volume and scale it to modern total ocean volume
        print('volumes scaled to modern ocean volume')
        box_vol *= MODERN_VOLUME/box_vol.sum()

    # atmosphere box
    box_vol[-1] = ATM_VOL

    # Conversion to GEOCLIM units
    box_vol *= VOLUME_CONVERSION_FACTOR

    # Write in file
    np.savetxt(outdir+vol_outfile, box_vol, fmt='%.12e', delimiter='\n')





###################################################################################################


## ----------------- ##
## Example for FOAM: ##
## ----------------- ##
#write_oceanic_input(co2=[560, 840, 1120, 1210, 1305, 1400, 1680, 1960, 2240, 4480, 8960],
#                    input_files=['560ppm.nc', '840ppm.nc', '1120ppm.nc', '1210ppm.nc', '1305ppm.nc',
#                                 '1400ppm.nc', '1680ppm.nc', '1960ppm.nc', '2240ppm.nc', '4480ppm.nc', '8960ppm.nc'],
#                    root='/home/piermafrost/Downloads/150Ma/merg150_VegDef_AdjCSol_EccN_ocean_',
#                    latitude='lat', z='lev', temperature='TEMP',
#                    y_weight='lw', z_weight='thickness')

## ------------------ ##
## Examples for IPSL: ##
## ------------------ ##
#
#write_oceanic_input(co2=[284.7],
#                    input_files=['CTRL-CM5A2_ANNCLIMO_oce.nc'], grid_file='coordinates.nc',
#                    root='./',
#                    latitude='nav_lat', z='deptht', temperature='thetao',
#                    horiz_weight='area', z_bounds='deptht_bounds')
#
#write_oceanic_input(co2=[1138.8],
#                    input_files=['4X/CPL-90Ma-ORB7a-TopoCorr_SE_8555_8654_1Y_grid_T.nc'],
#                    grid_file='PALEORCA2.90MaCorrected_grid.nc',
#                    temp_outfile='oceanic_temp_IPSL_90Ma_4PAL.dat',
#                    vol_outfile='oce_vol_IPSL-90Ma.dat',
#                    surf_outfile='oce_surf_IPSL-90Ma.dat',
#                    sedsurf_outfile='surf_sedi_IPSL-90Ma.dat',
#                    root='./',
#                    latitude='nav_lat', z='deptht', temperature='thetao',
#                    x_weight='e1t', y_weight='e2t', z_bounds='deptht_bounds')



# ========== #


# IPSL CTRL:
write_oceanic_input(co2=[284.7],
                    input_files=['piControl_SE_2750_2849_1Y_grid_T.nc'],
                    input_u_files=['piControl_SE_2750_2849_1Y_grid_U.nc'],
                    input_v_files=['piControl_SE_2750_2849_1Y_grid_V.nc'],
                    input_w_files=['piControl_SE_2750_2849_1Y_grid_W.nc'],
                    grid_file='coordinates.nc',
                    root='./',
                    outdir='IPSL-CM5A2_CTRL/',
                    latitude='nav_lat', z='deptht', temperature='thetao',
                    #u_flx='uocetr_eff', v_flx='vocetr_eff', w='wo', inverted_w=True,
                    u='uo', v='vo', w='wo', inverted_w=True,
                    horiz_weight='cell_area',
                    x_weight='e1t', y_weight='e2t', z_weight='e3t',
                    special_wrap='ORCA')


# IPSL 90Ma -- single climate:
#write_oceanic_input(co2=[4*284.7],
#                    input_files=['90Ma_oce_out/4X/CPL-90Ma-ORB7a-TopoCorr_SE_8555_8654_1Y_grid_T.nc'],
#                    grid_file='PALEORCA2.90MaCorrected_grid.nc',
#                    temp_outfile='oceanic_temp_90Ma_laugie_TopoCorr_4X_obl24.6_ecc0.015_pre0.dat',
#                    vol_outfile='oce_vol_IPSL-90Ma.dat',
#                    surf_outfile='oce_surf_IPSL-90Ma.dat',
#                    sedsurf_outfile='surf_sedi_IPSL-90Ma.dat',
#                    root='./',
#                    latitude='nav_lat', z='deptht', temperature='thetao',
#                    horiz_weight='srft', z_bounds='deptht_bnds')

