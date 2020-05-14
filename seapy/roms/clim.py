#!/usr/bin/env python
"""
  clim.py

  ROMS climatology utilities

  Written by Brian Powell on 08/15/15
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


import seapy
import numpy as np
import netCDF4


def gen_bry_clim(clim_file, grid, bry, clobber=False, cdl=None):
    """
    Taking the results of gen_ncks and interpolation, stitch together
    climatology files that were interpolated using only the boundary regions
    into a single climatology (with no data where interpolation wasn't
    performed).

    Parameters
    ----------
    clim_file: str,
        The name of the output climate file
    grid: seapy.model.grid or str,
        The output ROMS grid
    bry: dict,
        A dictionary prescribing the climatology file interpolated for each
        boundary side.
        {"west":filename, "south":filename}, ...}
    clobber: bool, optional
        If True, clobber any existing files and recreate. If False, use
        the existing file definition
    cdl: string, optional,
        Use the specified CDL file as the definition for the new
        netCDF file.

    Returns
    -------
        None
    """
    grid = seapy.model.asgrid(grid)

    # Grab the first dictionary record and use it to determine the number
    # of times in the new climatology file
    nc = netCDF4.Dataset(bry[list(bry.keys())[0]])
    reftime, time = seapy.roms.get_reftime(nc)
    times = nc.variables[time][:]
    nc.close()

    # Create the new climatology file
    ncout = seapy.roms.ncgen.create_clim(clim_file,
                                         eta_rho=grid.ln, xi_rho=grid.lm, s_rho=grid.n,
                                         reftime=reftime, clobber=clobber, cdl=cdl,
                                         title="stitched from boundary interpolation")
    ncout.variables["clim_time"][:] = times

    for side in bry:
        if bry[side] is None:
            continue
        ncin = netCDF4.Dataset(bry[side])
        for fld in seapy.roms.fields:
            idx = [np.s_[:] for i in range(seapy.roms.fields[fld]["dims"] + 1)]
            dat = ncin.variables[fld][:]
            shp = dat.shape
            if side == "west":
                idx[-1] = np.s_[:shp[-1]]
                pass
            elif side == "east":
                idx[-1] = np.s_[-shp[-1]:]
                pass
            elif side == "north":
                idx[-2] = np.s_[-shp[-2]:]
                pass
            elif side == "south":
                idx[-2] = np.s_[:shp[-2]]
                pass
            ncout.variables[fld][idx] = dat
            ncout.sync()
        ncin.close()
    ncout.close()
