#!/usr/bin/env python
"""
  analysis.py

  Methods to assist in the analysis of ROMS fields

  Written by Brian Powell on 05/24/15
  Copyright (c)2010--2025 University of Hawaii under the MIT-License.
"""

import numpy as np
from joblib import Parallel, delayed
import seapy
from rich.progress import track
from warnings import warn


def __find_surface_thread(grid, field, value, zeta, const_depth=False,
                          k_values=False, u_grid=False, v_grid=False):
    """
    Internal function to find a value in field_a and return the
    values from field_b at the same positions.
    """
    if zeta[0][0] is None:
        depth = grid.depth_rho
        if u_grid:
            depth = getattr(grid, 'depth_u', depth)
        elif v_grid:
            depth = getattr(grid, 'depth_v', depth)
    else:
        depth = seapy.roms.depth(grid.vtransform, grid.h, grid.hc,
                                 grid.s_rho, grid.cs_r, zeta)
        if u_grid:
            depth = seapy.model.rho2u(depth)
        elif v_grid:
            depth = seapy.model.rho2v(depth)

    # Set the variables based on what we are finding
    if const_depth:
        field_a, field_b = depth, field
    else:
        field_a, field_b = field, depth

    # Determine the upper and lower bounds of the value in the field
    tmp = np.ma.masked_equal(
        np.diff(((field_a - value) < 0).astype(np.intp), axis=0), 0)
    factor = -np.sign(np.mean(np.diff(field, axis=0))).astype(np.intp)

    # Determine the points of the upper bound and the lower bound
    bad = np.sum(tmp, axis=0).astype(bool)
    k_ones = np.arange(grid.n, dtype=np.intp)
    upper = (k_ones[:, np.newaxis, np.newaxis] ==
             np.argmax(np.abs(tmp), axis=0) + 1) * bad
    k_ones = np.arange(grid.n, dtype=np.intp) - factor
    lower = (k_ones[:, np.newaxis, np.newaxis] ==
             np.argmax(np.abs(tmp), axis=0) + 1) * bad

    # Now that we have the bounds, we can linearly interpolate to
    # find where the value lies
    u_a = np.sum(field_a * upper, axis=0)
    d_a = u_a - np.sum(field_a * lower, axis=0)
    d_z = (u_a - value) / d_a
    if k_values:
        return np.argmax(upper, axis=0) + factor * d_z

    # Calculate the values from field_b
    u_b = np.sum(field_b * upper, axis=0)
    d_b = u_b - np.sum(field_b * lower, axis=0)
    return u_b - d_b * d_z


def constant_depth(field, grid, depth, zeta=None, threads=2):
    """
    Find the values of a 3-D field at a constant depth for all times given.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to interpolate onto a constant depth level. If 4-D, it
        will calculate through time.
    grid : seapy.model.grid or string or list,
        Grid that defines the depths and stretching for the field given
    depth : float,
        Depth (in meters) to find all values
    zeta : ndarray, optional,
        ROMS zeta field corresponding to field if you wish to apply the SSH
        correction to the depth calculations.
    threads : int, optional,
        Number of threads to use for processing

    Returns
    -------
    nfield : ndarray,
        Values from ROMS field on the given constant depth
    """
    grid = seapy.model.asgrid(grid)
    field = np.ma.masked_invalid(field, copy=False)
    depth = depth if depth < 0 else -depth
    if depth is None or grid.depth_rho.min() > depth > grid.depth_rho.max():
        warn("Error: {:f} is out of range for the depth.".format(value))
        return
    if np.ndim(field) == 3:
        field = seapy.adddim(field)
    nt = field.shape[0]
    threads = np.minimum(nt, threads)
    if zeta is None:
        zeta = np.array(None * nt)[:, np.newaxis]
    if np.ndim(zeta) == 2:
        zeta = seapy.adddim(zeta, nt)

    v_grid = u_grid = False
    if field.shape[-2:] == grid.mask_u.shape:
        u_grid = True
    elif field.shape[-2:] == grid.mask_v.shape:
        v_grid = True

    return np.ma.array(Parallel(n_jobs=threads, verbose=0)
                       (delayed(__find_surface_thread)
                        (grid, field[i, :], depth, zeta[i, ...],
                         const_depth=True, u_grid=u_grid, v_grid=v_grid)
                        for i in range(nt)), copy=False)


def constant_value(field, grid, value, zeta=None, threads=2):
    """
    Find the depth of the value across the field. For example, find the depth
    of a given isopycnal if the field is density.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to interpolate onto a constant depth level. If 4-D, it
        will calculate through time.
    grid : seapy.model.grid or string or list,
        Grid that defines the depths and stretching for the field given
    value : float,
        Value to find the depths for in same units as the 'field'
    zeta : ndarray, optional,
        ROMS zeta field corresponding to field if you wish to apply the SSH
        correction to the depth calculations.
    threads : int, optional,
        Number of threads to use for processing

    Returns
    -------
    nfield : ndarray,
        Depths from ROMS field on the given value
    """
    grid = seapy.model.asgrid(grid)
    field = np.ma.masked_invalid(field, copy=False)
    if value is None or field.min() > value > field.max():
        warn("Error: {:f} is out of range for the field.".format(value))
        return
    if np.ndim(field) == 3:
        field = seapy.adddim(field)
    nt = field.shape[0]
    threads = np.minimum(nt, threads)
    if zeta is None:
        zeta = np.array([None] * nt)[:, np.newaxis]
    if np.ndim(zeta) == 2:
        zeta = seapy.adddim(zeta, nt)

    v_grid = u_grid = False
    if field.shape[-2:] == grid.mask_u.shape:
        u_grid = True
    elif field.shape[-2:] == grid.mask_v.shape:
        v_grid = True

    return np.ma.array(Parallel(n_jobs=threads, verbose=0)
                       (delayed(__find_surface_thread)
                        (grid, field[i, :], value, zeta[i, ...],
                         u_grid=u_grid, v_grid=v_grid)
                        for i in range(nt)), copy=False)


def constant_value_k(field, grid, value, zeta=None, threads=2):
    """
    Find the layer number of the value across the field. For example, find
    the k of a given isopycnal if the field is density.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to interpolate onto a constant depth level. If 4-D, it
        will calculate through time.
    grid : seapy.model.grid or string or list,
        Grid that defines the depths and stretching for the field given
    value : float,
        Value to find the depths for in same units as the 'field'
    zeta : ndarray, optional,
        ROMS zeta field corresponding to field if you wish to apply the SSH
        correction to the depth calculations.
    threads : int, optional,
        Number of threads to use for processing

    Returns
    -------
    nfield : ndarray,
        Depths from ROMS field on the given value
    """
    grid = seapy.model.asgrid(grid)
    field = np.ma.masked_invalid(field, copy=False)
    if value is None or field.min() > value > field.max():
        warn("Error: {:f} is out of range for the field.".format(value))
        return
    if np.ndim(field) == 3:
        field = seapy.adddim(field)
    nt = field.shape[0]
    threads = np.minimum(nt, threads)
    if zeta is None:
        zeta = np.zeros((nt, 1, 1))
    if np.ndim(zeta) == 2:
        zeta = seapy.adddim(zeta, nt)

    v_grid = u_grid = False
    if field.shape[-2:] == grid.mask_u:
        u_grid = True
    elif field.shape[-2:] == grid.mask_v:
        v_grid = True

    return np.ma.array(Parallel(n_jobs=threads, verbose=0)
                       (delayed(__find_surface_thread)
                        (grid, field[i, :], value, zeta[i, ...],
                         k_values=True, u_grid=u_grid, v_grid=v_grid)
                        for i in range(nt)), copy=False)


def gen_k_mask(N, kbot, ktop=None):
    """
    Given a k matrix, such as from constant_value_k(), generate
    a full-grid matrix that contains the weights of all values included.
    This allows you to perform multiplications to get proper contributions
    from the fractional layers.

    Parameters
    ----------
    N : int
        Number of layers present in the grid
    kbot : ndarray,
        A field of 'k' values that specify the bottom fractional layer.
    ktop : ndarray, optional
        A field of 'k' values that specify the top fractional layer. If not
        specified, the surface is assumed

    Returns
    -------
    nfield : ndarray,
        3-D field in size of k that is a weighting matrix. Value greater
        than zero is the weight of that cell. Value of zero means that cell
        is not used.

    """
    kbot = np.asarray(kbot)
    if ktop is None:
        ktop = N
    ktop = np.asarray(ktop)
    k_ones = np.arange(N, dtype=intp)
    dbot, _ = np.modf(kbot)
    dtop, _ = np.modf(ktop)
    fld = np.logical_and(
        k_ones[:, np.newaxis, np.newaxis] >= kbot.astype(int),
        k_ones[:, np.newaxis, np.newaxis] < np.ceil(ktop).astype(int)).astype(float)
    bfrac = (k_ones[:, np.newaxis, np.newaxis] == kbot.astype(int)).astype(float)
    tfrac = (k_ones[:, np.newaxis, np.newaxis] == ktop.astype(int)).astype(float)

    return fld - bfrac * dbot - tfrac * (1 - dtop)


def depth_sum(field, thickness, bottom, top=0, partial=False, average=False):
    """
    Compute the depth-integrated field between the specified bottom and top
    depths. Because a grid cell represents a volume, the thickness is used
    rather than depth. This provides the most accurate integration.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to integrate from a depth level. Must be
        three-dimensional array (single time).
    thickness : ndarray,
        3-D array of the thickness of each grid cell in field
    bottom : float,
        Depth (in meters) to integrate from
    top : float, [optional]
        Depth (in meters) to integrate to defaults to surface
    partial: boolean, [optional]
        If True, produce results from grid locations that don't cover the
        full range of depths specified (i.e., bottom is 30, but the call
        was to integrate between 50 and 20m). Default is False.
    average: boolean, [optional]
        When False [default], depth integrate the field. If True, calculate
        the depth average.

    Returns
    -------
    ndarray,
        Values from depth integrated ROMS field
    """
    bottom = bottom if bottom > 0 else -bottom
    top = top if top > 0 else -top
    if bottom < top:
        bottom, top = top, bottom

    # Check that the span and field are consistent
    idim = field.ndim - thickness.ndim
    if field.shape[idim:] != thickness.shape:
        raise ValueError("field and thickness cannot be broadcast together: " +
                         f"{field.shape} vs {thickness.shape}")

    # Calculate the depth span of each cell
    span = np.cumsum(thickness[::-1, ...], axis=0)[::-1, ...]

    # Set up arrays for indexing
    surf = span.shape[0]
    k_ones = np.arange(span.shape[0], dtype=intp)
    sl1 = np.arange(span.shape[1], dtype=intp)[:, None]
    sl2 = np.arange(span.shape[2], dtype=intp)[None, :]

    # Helper function to find fractional minimums
    def _find_fraction_layer0(f):
        # Identify problem locations, where there is no zero crossing
        shallow = f.min(axis=0) > 0
        deep = f.max(axis=0) <= 0
        mask = np.argmin(np.ma.masked_equal(k_ones[:, np.newaxis, np.newaxis] *
                                            (f < 0), 0), axis=0)
        mask1 = np.maximum(0, mask - 1)
        return mask - 1 + f[mask1, sl1, sl2] / thickness[mask1, sl1, sl2], \
            shallow, deep

    # Find the bottom k levels
    kbot, shbot, dpbot = _find_fraction_layer0(span - bottom)
    # If the bottom is too shallow, we can just use the surface
    kbot[shbot] = surf - 1
    # If the bottom is too deep, use the bottom
    kbot[dpbot] = 0

    # Find the top k levels
    ktop = np.array([surf])
    dptop = np.full(kbot.shape, False)
    if top != 0:
        ktop, shtop, dptop = _find_fraction_layer0(span - top)
        # If the top is too shallow, use the surface
        ktop[shtop] = surf

    # Build the mask
    mask = gen_k_mask(surf, kbot, ktop)

    # Set mask to nan for locations with a water depth above the top
    mask[:, dptop] = np.nan
    # Set mask to nan for locations with a water depth that only partially
    # covers the specified region unless the user requests partial points
    if not partial:
        mask[:, dpbot] = np.nan

    # Do the integration
    if average:
        return np.ma.masked_invalid(np.sum(field * mask * thickness, axis=idim) /
                                    np.sum(thickness * mask, axis=0), copy=False)
    else:
        return np.ma.masked_invalid(np.sum(field * mask * thickness, axis=idim))


def depth_average(field, thickness, bottom, top=0, partial=False):
    """
    Compute the depth-averaged field between the specified bottom and top
    depths. Because a grid cell represents a volume, the thickness is used
    rather than depth. This provides the most accurate integration.

    Parameters
    ----------
    field : ndarray,
        ROMS 3-D field to integrate from a depth level. Must be
        three-dimensional array (single time).
    thickness : ndarray,
        3-D array of the thickness of each grid cell in field
    bottom : float,
        Depth (in meters) to integrate from
    top : float, [optional]
        Depth (in meters) to integrate to defaults to surface
    partial: boolean, [optional]
        If True, produce results from grid locations that don't cover the
        full range of depths specified (i.e., bottom is 30, but the call
        was to integrate between 50 and 20m). Default is False.

    Returns
    -------
    ndarray,
        Values from depth integrated ROMS field
    """
    return depth_sum(field, thickness, bottom, top, partial, average=True)


def transect(lon, lat, depth, data, nx=200, nz=40, z=None):
    """
    Generate an equidistant transect from data at varying depths. Can be
    used to plot a slice of model or observational data.

    Parameters
    ----------
    lat: array
        n-dimensional array with latitude of points
    lon: array
        n-dimensional array with longitude of points
    depth: array
        [k,n] dimensional array of the depths for all k-layers at each n point
    data: array
        [k,n] dimensional array of values
    nx: int, optional
        number of horizontal points desired in the transect
    nz: int, optional
        number of vertical points desired in the transect
    z: array, optional
        list of depths to use if you do not want equidistant depths

    Returns
    -------
    x: array
        x-location values in [m] along transect
    z: array
        depth values in [m] of the new transect
    vals: np.ma.array
        data values of the new transect with masked values

    Examples
    --------
    Generate series of transects from ROMS output

    >>> nc = seapy.netcdf('roms_his.nc')
    >>> grid = seapy.model.asgrid(nc)
    >>> data = nc.variables['salt'][:,:,150,:]
    >>> shp = (data.shape[0], 50, 400)
    >>> transect = np.zeros(shp)
    >>> for i in range(shp[0]):
    >>>     x, z, d = \
    >>>          seapy.roms.analysis.transect(grid.lon_rho[150,:],
    >>>                                       grid.lat_rho[150,:],
    >>>                                       grid.depth_rho[:,150,:],
    >>>                                       data[i,:,:], nx=shp[2],
    >>>                                       nz=shp[1])
    >>>     transect[i,:,:] = d.filled(np.nan)
    >>> nc.close()
    >>> plt.pcolormesh(x/1000, z, transect[0, :, :])
    """
    from scipy.interpolate import griddata
    depth = np.atleast_2d(depth)
    data = np.ma.atleast_2d(data).filled(np.mean(data))
    lon = np.atleast_1d(lon)
    lat = np.atleast_1d(lat)

    # Generate the depths
    depth[depth > 0] *= -1
    if z is None:
        z = np.linspace(depth.min() - 2, depth.max(), nz)
    else:
        z[z > 0] *= -1
        nz = len(z)
    dz = np.abs(np.diff(z).mean())

    # Determine the distance between points and the weighting to apply
    dist = np.hstack(([0], seapy.earth_distance(
        lon[0], lat[0], lon[1:], lat[1:])))
    dx = np.diff(dist).mean()
    zscale = np.maximum(1, 10**int(np.log10(dx / dz)))
    dx /= zscale
    x = np.linspace(0, dist.max(), nx)

    # All arrays have to be the same size
    xx, zz = np.meshgrid(x / zscale, z)

    # For the source data, we don't want to extrpolate,
    # so make the data go from the surface to twice its
    # depth.
    zl = np.argsort(depth[:, 0])
    dep = np.vstack((np.ones((1, depth.shape[1])) * 2 * depth.min(),
                     depth[zl, :],
                     np.zeros((1, depth.shape[1]))))

    # repeat the same data at the top and bottom
    dat = np.vstack((data[zl[0], :], data[zl],
                     data[zl[-1], :]))
    dist = np.tile(dist.T, [dep.shape[0], 1]) / zscale

    # Find the bottom indices to create a mask for nodata/land
    idx = np.interp(xx[0, :], dist[0, :],
                    np.interp(depth.min(axis=0), z,
                              np.arange(nz))).astype(int)
    mask = np.arange(nz)[:, np.newaxis] <= idx

    # Interpolate
    ndat = np.ma.array(griddata(
        (dist.ravel(), dep.ravel()), dat.ravel(), (xx.ravel(), zz.ravel()),
        method='cubic').reshape(xx.shape), mask=mask)

    # Return everything
    return x, z, ndat


def gen_std_i(roms_file, std_file, std_window=5, pad=1, skip=30, fields=None):
    """
    Create a std file for the given ocean fields. This std file can be used
    for initial conditions constraint in 4D-Var. This requires a long-term
    model spinup file from which to compute the standard deviation.

    Parameters
    ----------
    roms_file: string or list of strings,
        The ROMS (history or average) file from which to compute the std. If
        it is a list of strings, a netCDF4.MFDataset is opened instead.
    std_file: string,
        The name of the file to store the standard deviations fields
    std_window: int,
        The size of the window (in number of records) to compute the std over
    pad: int,
        How much to pad each side of the window for overlap. For example,
        std_window=10 and pad=2 would give a total window of 14 with 2 records
        used in the prior window and 2 in the post window as well.
    skip: int,
        How many records to skip at the beginning of the file
    fields: list of str,
        The fields to compute std for. Default is to use the ROMS prognostic
        variables.

    Returns
    -------
        None
    """
    # Create the fields to process
    if fields is None:
        fields = set(seapy.roms.fields)

    # Open the ROMS info
    grid = seapy.model.asgrid(roms_file)
    nc = seapy.netcdf(roms_file)

    # Filter the fields for the ones in the ROMS file
    fields = set(nc.variables).intersection(fields)

    # Build the output file
    epoch, time_var = seapy.roms.get_reftime(nc)
    time = nc.variables[time_var][:]
    ncout = seapy.roms.ncgen.create_da_ini_std(std_file,
                                               eta_rho=grid.ln, xi_rho=grid.lm,
                                               s_rho=grid.n,
                                               reftime=epoch,
                                               title="std from " + str(roms_file))
    grid.to_netcdf(ncout)

    # If there are any fields that are not in the standard output file,
    # add them to the output file
    for f in fields.difference(ncout.variables):
        ncout.createVariable(f, np.float32,
                             ('ocean_time', "s_rho", "eta_rho", "xi_rho"))

    # Loop over the time with the variance window:
    time_list = np.arange(skip + pad, len(time) - std_window - pad, std_window)
    for n, t in track(enumerate(time_list), total=len(time_list),
                      description="evaluate time window"):
        idx = np.arange(t - pad, t + std_window + pad)
        ncout.variables[time_var][n] = np.mean(time[idx])
        for v in fields:
            dat = nc.variables[v][idx, :].std(axis=0)
            dat[dat > 10] = 0.0
            ncout.variables[v][n, :] = dat
        ncout.sync()
    ncout.close()
    nc.close()


def gen_std_f(roms_file, std_file, records=None, fields=None):
    """
    Create a std file for the given atmospheric forcing fields. This std
    file can be used for the forcing constraint in 4D-Var. This requires a
    long-term model spinup file from which to compute the standard deviation.

    Parameters
    ----------
    roms_file: string or list of strings,
        The ROMS (history or average) file from which to compute the std. If
        it is a list of strings, a netCDF4.MFDataset is opened instead.
    std_file: string,
        The name of the file to store the standard deviations fields
    records: ndarray,
        List of records to perform the std over. These records are used to
        avoid the solar diurnal cycles in the fields.
    fields: list of str,
        The fields to compute std for. Default is to use the ROMS atmospheric
        variables (sustr, svstr, shflux, ssflux).

    Returns
    -------
        None
    """
    # Create the fields to process
    if fields is None:
        fields = set(["sustr", "svstr", "shflux", "ssflux"])

    # Open the ROMS info
    grid = seapy.model.asgrid(roms_file)
    nc = seapy.netcdf(roms_file)

    # Filter the fields for the ones in the ROMS file
    fields = set(nc.variables).intersection(fields)

    # Build the output file
    epoch, time_var = seapy.roms.get_reftime(nc)
    time = nc.variables[time_var][:]
    ncout = seapy.roms.ncgen.create_da_frc_std(std_file,
                                               eta_rho=grid.ln, xi_rho=grid.lm, s_rho=grid.n,
                                               reftime=epoch, title="std from " + str(roms_file))
    grid.to_netcdf(ncout)

    # Set the records
    if records is None:
        records = np.arange(len(time))
    else:
        records = np.atleast_1d(records)
        records = records[records <= len(time)]

    # If there are any fields that are not part of the standard, add them
    # to the output file
    for f in fields.difference(ncout.variables):
        ncout.createVariable(f, np.float32,
                             ('ocean_time', "eta_rho", "xi_rho"))

    # Loop over the time with the variance window:
    ncout.variables[time_var][:] = np.mean(time[records])
    for v in fields:
        dat = nc.variables[v][records, :].std(axis=0)
        ncout.variables[v][0, :] = dat
        ncout.sync()
    ncout.close()
    nc.close()


def plot_obs_spatial(obs, type='zeta', prov=None, time=None, depth=0,
                     gridcoord=False, error=False, **kwargs):
    """
    Create a surface plot of the observations.

    Parameters
    ----------
    obs: filename, list, or observation class
        The observations to use for plotting
    type: string or int,
        The type of observation to plot ('zeta', 'temp', 'salt', etc.)
    prov: string or int,
        The provenance of the observations to plot
    time: ndarray,
        The times of the observations to plot
    depth: float,
        The depth of the obs to plot over the spatial region
    gridcoord: bool,
        If True, plot on grid coordinates. If False [default] plot on lat/lon
    error: bool,
        If True plot the errors rather than values. Default is False.
    **kwargs: keywords
        Passed to matplotlib.pyplot.scatter

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt

    obs = seapy.roms.obs.asobs(obs)
    otype = seapy.roms.obs.astype(type)
    if prov is not None:
        prov = seapy.roms.obs.asprovenance(prov)
    if time is not None:
        time = np.atleast_1d(time)

    # Search the obs for the user
    if prov is not None:
        idx = np.where(np.logical_and.reduce((
            obs.type == otype,
            obs.provenance == prov,
            np.logical_or(obs.z == 0, obs.depth == depth))))[0]

    else:
        idx = np.where(np.logical_and(
            obs.type == otype,
            np.logical_or(obs.z == 0, obs.depth == depth)))[0]

    # If there is a time specific condition, find the sets
    if time is not None:
        idx = idx[np.in1d(obs.time[idx], time)]

    # If we don't have anything to plot, return
    if not idx.any():
        return

    # Plot it up
    if not kwargs:
        kwargs = {'s': 30, 'alpha': 0.8, 'linewidths': (0, 0)}
    if gridcoord:
        x = obs.x
        y = obs.y
    else:
        x = obs.lon
        y = obs.lat
    val = obs.value if not error else np.sqrt(obs.error)
    plt.scatter(x[idx], y[idx], c=val[idx], **kwargs)
    plt.colorbar()


def plot_obs_profile(obs, type='temp', prov=None, time=None,
                     gridcoord=False, error=False, **kwargs):
    """
    Create a sub-surface profile plot of the observations.

    Parameters
    ----------
    obs: filename, list, or observation class
        The observations to use for plotting
    type: string or int,
        The type of observation to plot ('zeta', 'temp', 'salt', etc.)
    prov: string or int,
        The provenance of the observations to plot
    time: ndarray,
        The times of the observations to plot
    gridcoord: bool,
        If True, plot on grid coordinates. If False [default] plot on lat/lon
    error: bool,
        If True plot the errors rather than values. Default is False.
    **kwargs: keywords
        Passed to matplotlib.pyplot.scatter

    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt

    obs = seapy.roms.obs.asobs(obs)
    otype = seapy.roms.obs.astype(type)
    if prov is not None:
        prov = seapy.roms.obs.asprovenance(prov)
    if time is not None:
        time = np.atleast_1d(time)

    # Search the obs for the user
    if prov is not None:
        idx = np.where(np.logical_and.reduce((
            obs.type == otype,
            obs.provenance == prov,
            np.logical_or(obs.z < 0, obs.depth < 0))))[0]

    else:
        idx = np.where(np.logical_and(
            obs.type == otype,
            np.logical_or(obs.z < 0, obs.depth < 0)))[0]

    # If there is a time specific condition, find the sets
    if time is not None:
        idx = idx[np.in1d(obs.time[idx], time)]

    # If we don't have anything to plot, return
    if not idx.any():
        return

    # Plot it up
    if gridcoord:
        dep = obs.z if np.mean(obs.z[idx] > 0) else obs.depth
    else:
        dep = obs.z if np.mean(obs.z[idx] < 0) else obs.depth
    val = obs.value if not error else np.sqrt(obs.error)
    plt.plot(val[idx], dep[idx], 'k+', **kwargs)
