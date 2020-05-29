#!/usr/bin/env python
"""
  grid

  This module handles general model grid information, whether from ROMS or
  other models; however, it is mostly geared towards ROMS

  Written by Brian Powell on 10/09/13
  Copyright (c)2020 University of Hawaii under the MIT-License.

  **Examples**

  >>> grid = seapy.model.asgrid("grid_file.nc")

"""


import os
import re
import seapy
import numpy as np
import scipy.spatial
import matplotlib.path
import netCDF4
from warnings import warn


def asgrid(grid):
    """
    Return either an existing or new grid object. This decorator will ensure that
    the variable being used is a seapy.model.grid. If it is not, it will attempt
    to construct a new grid with the variable passed.

    Parameters
    ----------
    grid: string, list, netCDF4 Dataset, or model.seapy.grid
        Input variable to cast. If it is already a grid, it will return it;
        otherwise, it attempts to construct a new grid.

    Returns
    -------
    seapy.model.grid

    """
    if grid is None:
        raise AttributeError("No grid was specified")
    if isinstance(grid, seapy.model.grid):
        return grid
    if isinstance(grid, netCDF4._netCDF4.Dataset):
        return seapy.model.grid(nc=grid)
    else:
        return seapy.model.grid(filename=grid)


class grid:

    def __init__(self, filename=None, nc=None, lat=None, lon=None, z=None,
                 depths=True, cgrid=False):
        """
            Class to wrap around a numerical model grid for oceanography.
            It attempts to track latitude, longitude, z, and other
            parameters. A grid can be constructed by specifying a filename or
            by specifying lat, lon, and z.

            Parameters
            ----------
            filename: filename or list, optional
              name to load to build data structure
                or
            nc: netCDF4 Dataset, optional
              If a file is already open, pass the reference.
            lat: ndarray,
                latitude values of grid
            lon: ndarray,
                longitude values of grid
            z : ndarray,
                z-level depths of grid

            Options
            -------
            depths: ndarray,
                Set the depths of the grid [True]
            cgrid: bool,
                Whether the grid is an Arakawa C-Grid [False]
        """
        self.filename = filename
        self.cgrid = cgrid
        self._nc = nc

        if (self.filename or self._nc) is not None:
            self._initfile()
            self._isroms = True if \
                (len(list(set(("s_rho", "pm", "pn", "theta_s", "theta_b",
                               "vtransform", "vstretching")).intersection(
                    set(self.__dict__)))) > 0) else False
            self.cgrid = True if self._isroms else self.cgrid
        else:
            self._nc = None
            self._isroms = False
            self.lat_rho = lat
            self.lon_rho = lon
            self.z = z
            self.cgrid = False
        self._verify_shape()
        if depths:
            self.set_dims()
            self.set_depth()
            self.set_thickness()
            self.set_mask_h()
        self.ijinterp = None
        self.llinterp = None

    def _initfile(self):
        """
        Using an input file, try to load as much information
        as can be found in the given file.

        Parameters
        ----------
        None

        Returns
        -------
        None : sets attributes in grid

        """
        # Define a dictionary to go through and convert netcdf variables
        # to internal class attributes
        gvars = {"lat_rho": ["lat_rho", "lat", "latitude", "y_rho", "geolat_t"],
                 "lon_rho": ["lon_rho", "lon", "longitude", "x_rho", "geolon_t"],
                 "lat_u": ["lat_u", "y_u", "geolat_u"],
                 "lon_u": ["lon_u", "x_u", "geolon_u"],
                 "lat_v": ["lat_v", "y_v", "geolat_u"],
                 "lon_v": ["lon_v", "x_v", "geolon_u"],
                 "mask_rho": ["mask_rho", "mask"],
                 "mask_u": ["mask_u"],
                 "mask_v": ["mask_v"],
                 "angle": ["angle"],
                 "h": ["h"],
                 "n": ["n"],
                 "theta_s": ["theta_s"],
                 "theta_b": ["theta_b"],
                 "tcline": ["tcline"],
                 "hc": ["hc"],
                 "vtransform": ["vtransform"],
                 "vstretching": ["vstretching"],
                 "s_rho": ["s_rho"],
                 "cs_r": ["cs_r"],
                 "f": ["f"],
                 "pm": ["pm"],
                 "pn": ["pn"],
                 "z": ["z", "depth", "lev", "st_ocean"],
                 "wtype_grid": ["mask_rho"],
                 "rdrag": ["rdrag"],
                 "rdrag2": ["rdrag2"],
                 "diff_factor": ["diff_factor"],
                 "visc_factor": ["visc_factor"]
                 }

        # Open the file
        close = False
        if self._nc is None:
            close = True
            self._nc = seapy.netcdf(self.filename)
        try:
            self.name = re.search("[^\.]*",
                                  os.path.basename(self.filename)).group()
        except:
            self.name = "untitled"
        self.key = {}
        ncvars = {v.lower(): v for v in self._nc.variables.keys()}
        for var in gvars:
            for inp in gvars[var]:
                if inp in ncvars:
                    self.key[var] = inp
                    self.__dict__[var] = self._nc.variables[ncvars[inp]][:]
                    break

        if close:
            # Close the file
            self._nc.close()
            self._nc = None

    def _verify_shape(self):
        """
        Verify the dimensionality of the system, create variables that
        can be generated from the others if they aren't already loaded

        Parameters
        ----------
        None

        Returns
        -------
        None : sets attributes in grid
        """
        # Check that we have the minimum required data
        if ("lat_rho" or "lon_rho") not in self.__dict__:
            raise AttributeError(
                "grid does not have attribute lat_rho or lon_rho")

        # Check that it is formatted into 2-D
        self.spatial_dims = self.lat_rho.ndim
        if self.lat_rho.ndim == 1 and self.lon_rho.ndim == 1:
            [self.lon_rho, self.lat_rho] = np.meshgrid(self.lon_rho,
                                                       self.lat_rho)

        # Compute the dimensions
        self.ln = int(self.lat_rho.shape[0])
        self.lm = int(self.lat_rho.shape[1])
        self.shape = (self.ln, self.lm)
        if self.cgrid:
            self.shape_u = (self.ln, self.lm - 1)
            self.shape_v = (self.ln - 1, self.lm)
        else:
            self.shape_u = self.shape_v = self.shape

    def __repr__(self):
        return "{:s}: {:d}x{:d}x{:d}".format("C-Grid" if self.cgrid
                                             else "A-Grid", self.n, self.ln, self.lm)

    def __str__(self):
        return "\n".join((self.filename if self.filename else "Constructed",
                          "{:d}x{:d}x{:d}: {:s} with {:s}".format(
                              self.n, self.ln, self.lm,
                              "C-Grid" if self.cgrid else "A-Grid",
                              "S-level" if self._isroms else "Z-Level"),
                          "Available: " + ",".join(sorted(
                              list(self.__dict__.keys())))))

    def east(self):
        """
        Test the longitude convention of the grid. If there are negative
        values, then east is False. If there are only positive values then
        assume that east is True.

        Parameters
        ----------
        None

        Returns
        -------
        east : bool,
            True - The convention is all positive values to the east
            False - The convention is positive values east and negative west
        """
        return np.min(self.lon_rho > 0)

    def set_east(self, east=False):
        """
        When working with various other grids, we may want the longitudes
        to be consistent. This can be changed by setting the east to be
        either True or False. If False, then longitudes will be positive east
        and negative west. If True, only positive east.

        Parameters
        ----------
        east : bool,
            True - all longitudes are positive
            False - longitudes are positive east and negative west

        Returns
        -------
        None : sets attributes in grid
        """
        try:
            if east:
                self.lon_rho[self.lon_rho < 0] += 360.0
                self.lon_u[self.lon_u < 0] += 360.0
                self.lon_v[self.lon_v < 0] += 360.0
            else:
                self.lon_rho[self.lon_rho > 180] -= 360.0
                self.lon_u[self.lon_u > 180] -= 360.0
                self.lon_v[self.lon_v > 180] -= 360.0
        except:
            pass

    def set_dims(self):
        """
        Compute the dimension attributes of the grid based upon the information provided.

        Parameters
        ----------
        None

        Returns
        -------
        None : sets attributes in grid
        """
        # If C-Grid, set the dimensions for consistency
        if self.cgrid:
            self.eta_rho = self.ln
            self.eta_u = self.ln
            self.eta_v = self.ln - 1
            self.xi_rho = self.lm
            self.xi_u = self.lm - 1
            self.xi_v = self.lm

        # Set the number of layers
        if "n" not in self.__dict__:
            if "s_rho" in self.__dict__:
                self.n = int(self.s_rho.size)
            elif "z" in self.__dict__:
                self.n = int(self.z.size)
            else:
                self.n = 1
                self.z = np.zeros(self.lat_rho.shape)
        else:
            self.n = int(self.n)

        # Generate the u- and v-grids
        if ("lat_u" or "lon_u") not in self.__dict__:
            if self.cgrid:
                self.lat_u = 0.5 * \
                    (self.lat_rho[:, 1:] - self.lat_rho[:, 0:-1])
                self.lon_u = 0.5 * \
                    (self.lon_rho[:, 1:] - self.lon_rho[:, 0:-1])
            else:
                self.lat_u = self.lat_rho
                self.lon_u = self.lon_rho
        if ("lat_v" or "lon_v") not in self.__dict__:
            if self.cgrid:
                self.lat_v = 0.5 * \
                    (self.lat_rho[1:, :] - self.lat_rho[0:-1, :])
                self.lon_v = 0.5 * \
                    (self.lon_rho[1:, :] - self.lon_rho[0:-1, :])
            else:
                self.lat_v = self.lat_rho
                self.lon_v = self.lon_rho
        if "mask_rho" in self.__dict__:
            if "mask_u" not in self.__dict__:
                if self.cgrid:
                    self.mask_u = self.mask_rho[:, 1:] * self.mask_rho[:, 0:-1]
                else:
                    self.mask_u = self.mask_rho
            if "mask_v" not in self.__dict__:
                if self.cgrid:
                    self.mask_v = self.mask_rho[1:, :] * self.mask_rho[0:-1, :]
                else:
                    self.mask_v = self.mask_rho

        # Compute the resolution
        if "pm" in self.__dict__:
            self.dm = 1.0 / self.pm
        else:
            self.dm = np.ones(self.lon_rho.shape, dtype=np.float32)
            self.dm[:, 0:-1] = seapy.earth_distance(self.lon_rho[:, 1:],
                                                    self.lat_rho[:, 1:],
                                                    self.lon_rho[:, 0:-1],
                                                    self.lat_rho[:, 0:-1]).astype(np.float32)
            self.dm[:, -1] = self.dm[:, -2]
            self.pm = 1.0 / self.dm
        if "pn" in self.__dict__:
            self.dn = 1.0 / self.pn
        else:
            self.dn = np.ones(self.lat_rho.shape, dtype=np.float32)
            self.dn[0:-1, :] = seapy.earth_distance(self.lon_rho[1:, :],
                                                    self.lat_rho[1:, :],
                                                    self.lon_rho[0:-1, :],
                                                    self.lat_rho[0:-1, :]).astype(np.float32)
            self.dn[-1, :] = self.dn[-2, :]
            self.pn = 1.0 / self.dn

        # Compute the Coriolis
        if "f" not in self.__dict__:
            omega = 2 * np.pi * seapy.secs2day
            self.f = 2 * omega * np.sin(np.radians(self.lat_rho))

        # Set the grid index coordinates
        self.I, self.J = np.meshgrid(
            np.arange(0, self.lm), np.arange(0, self.ln))

    def set_mask_h(self, fld=None):
        """
        Compute the mask and h array from a z-level model

        Parameters
        ----------
        fld : np.array
            3D array of values (such as temperature) to analyze to determine
            where the bottom and land lie

        Returns
        -------
        None : sets mask and h attributes in grid

        """
        if hasattr(self, "mask_rho") or self.cgrid:
            return
        if fld is None and self.filename is not None:
            if self._nc is None:
                self._nc = seapy.netcdf(self.filename)

            # Try to load a field from the file
            for f in ["temp", "temperature", "water_temp", "fed"]:
                if f in self._nc.variables:
                    fld = self._nc.variables[f][0, :, :, :]
                    fld = np.ma.array(fld, mask=np.isnan(fld))
                    break

            # Close the file
            self._nc.close()

        # If we don't have a field to examine, then we cannot compute the
        # mask and bathymetry
        if fld is None:
            warn("Missing 3D field to evaluate.")
            return

        # Next, we go over the field to examine the depths and mask
        self.h = np.zeros(self.lat_rho.shape)
        self.mask_rho = np.zeros(self.lat_rho.shape)
        for k in range(self.z.size):
            water = np.nonzero(np.logical_not(fld.mask[k, :, :]))
            self.h[water] = self.z[k]
            if k == 0:
                self.mask_rho[water] = 1.0
        self.mask_u = self.mask_v = self.mask_rho

    def set_depth(self, force=False):
        """
        Compute the depth of each cell for the model grid.

        Parameters
        ----------
        force : boolean, default False
                If True, force the update of the depths

        Returns
        -------
        None : sets depth attributes in grid
        """
        try:
            if self._isroms:
                if "s_rho" not in self.__dict__ or \
                   "cs_r" not in self.__dict__ or force:
                    self.s_rho, self.cs_r = seapy.roms.stretching(
                        self.vstretching, self.theta_s, self.theta_b,
                        self.hc, self.n)
                self.depth_rho = seapy.roms.depth(
                    self.vtransform, self.h, self.hc, self.s_rho, self.cs_r)
                self.depth_u = seapy.model.rho2u(self.depth_rho).filled(0)
                self.depth_v = seapy.model.rho2v(self.depth_rho).filled(0)
            else:
                d = self.z.copy()
                l = np.nonzero(d > 0)
                d[l] = -d[l]
                if self.n > 1:
                    self.depth_rho = np.kron(np.kron(
                        d, np.ones(self.lon_rho.shape[1])),
                        np.ones(self.lon_rho.shape[0])).reshape(
                        [self.z.size, self.lon_rho.shape[0],
                         self.lon_rho.shape[1]])
                else:
                    self.depth_rho = self.z
                if self.cgrid:
                    self.depth_u = seapy.model.rho2u(self.depth_rho).filled(0)
                    self.depth_v = seapy.model.rho2v(self.depth_rho).filled(0)
                else:
                    self.depth_u = self.depth_rho
                    self.depth_v = self.depth_rho
        except (AttributeError, ValueError):
            warn("could not compute grid depths.")
            pass

    def set_thickness(self):
        """
        Compute the thickness of each cell for the model grid.

        Parameters
        ----------
        None

        Returns
        -------
        None : sets thick attributes in grid
        """
        if "n" not in self.__dict__:
            self.set_dims()
        if self.n == 1:
            return
        try:
            if self._isroms:
                s_w, cs_w = seapy.roms.stretching(
                    self.vstretching, self.theta_s, self.theta_b, self.hc,
                    self.n, w_grid=True)
                self.thick_rho = seapy.roms.thickness(
                    self.vtransform, self.h, self.hc, s_w, cs_w)
                self.thick_u = seapy.model.rho2u(self.thick_rho)
                self.thick_v = seapy.model.rho2v(self.thick_rho)
            else:
                d = np.abs(self.z.copy())
                w = d * 0
                # Check which way the depths are going
                if d[0] < d[-1]:
                    w[0] = d[0]
                    w[1:] = d[1:] - d[0:-1]
                else:
                    w[-1] = d[-1]
                    w[0:-1] = d[0:-1] - d[1:]

                self.thick_rho = np.kron(np.kron(w,
                                                 np.ones(self.lon_rho.shape[1])),
                                         np.ones(self.lon_rho.shape[0])).reshape(
                    [self.z.size, self.lon_rho.shape[0],
                     self.lon_rho.shape[1]])
                if self.cgrid:
                    self.thick_u = seapy.model.rho2u(self.thick_rho)
                    self.thick_v = seapy.model.rho2v(self.thick_rho)
                else:
                    self.thick_u = self.thick_rho
                    self.thick_v = self.thick_rho
        except AttributeError:
            warn("could not compute grid thicknesses.")
            pass

    def plot_trace(self, basemap=None, **kwargs):
        """
        Trace the boundary of the grid onto a map projection

        Parameters
        ----------
        basemap: basemap instance
            The basemap instance to use for drawing
        **kwargs: optional
            Arguments to pass to the plot routine

        Returns
        -------
        None
        """
        lon = np.concatenate([self.lon_rho[0, :], self.lon_rho[:, -1],
                              self.lon_rho[-1, ::-1], self.lon_rho[::-1, 0]])
        lat = np.concatenate([self.lat_rho[0, :], self.lat_rho[:, -1],
                              self.lat_rho[-1, ::-1], self.lat_rho[::-1, 0]])
        if basemap:
            x, y = basemap(lon, lat)
            basemap.plot(x, y, **kwargs)
        else:
            from matplotlib import pyplot
            pyplot.plot(lon, lat, **kwargs)

    def plot_depths(self, row=None, col=None, ax=None):
        """
        Plot the depths of a model grid along a row or column transect.
        If the bathymetry is known, it is plotted also.

        Parameters
        ----------
        row : int, optional
          The row number to plot
        col : int, optional
          The column number to plot
        ax : matplotlib.axes, optional
          The axes to use for the figure

        Returns
        -------
        ax : matplotlib.axes
          The axes containing the plot
        """
        import matplotlib.pyplot as plt

        # Create the axes if we don't have any
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # ax.set_bg_color('darkseagreen')

        # Get the data
        if row:
            sz = np.s_[:, row, :]
            s = np.s_[row, :]
            x = self.lon_rho[s]
            label = "Longitude"
        elif col:
            sz = np.s_[:, :, col]
            s = np.s_[:, col]
            x = self.lat_rho[s]
            label = "Latitude"
        else:
            warn("You must specify a row or column")
            return

        # If it is ROMS, we should plot the top and bottom of the cells
        if self._isroms:
            sr, csr = seapy.roms.stretching(
                self.vstretching, self.theta_s, self.theta_b,
                self.hc, self.n, w_grid=True)
            dep = np.ma.masked_where(seapy.adddim(self.mask_rho[s],
                                                  self.n + 1) == 0,
                                     seapy.roms.depth(self.vtransform,
                                                      self.h[s], self.hc,
                                                      sr, csr,
                                                      w_grid=True))
        else:
            dep = np.ma.masked_where(seapy.adddim(self.mask_rho[s],
                                                  self.n) == 0,
                                     self.depth_rho[sz])
        h = -self.h[s]

        # Begin with the bathymetric data
        ax.fill_between(x, h, np.min(h), facecolor="darkseagreen",
                        interpolate=True)

        # Plot the layers
        ax.plot(x, dep.T, color="grey")

        # Labels
        ax.set_xlabel(label + " [deg]")
        ax.set_ylabel("Depth [m]")

        # Make it tight
        plt.autoscale(ax, tight=True)

        return ax

    def to_netcdf(self, nc):
        """
        Write all available grid information into the records present in the
        netcdf file.  This is used to pre-fill boundary, initial, etc. files
        that require some of the grid information.

        Parameters
        ----------
        nc : netCDF4
            File to fill all known records from the grid information

        Returns
        -------
        None

        """
        for var in nc.variables:
            if hasattr(self, var.lower()):
                nc.variables[var][:] = getattr(self, var.lower())

    def nearest(self, lon, lat, grid="rho"):
        """
        Find the indices nearest to each point in the given list of
        longitudes and latitudes.

        Parameters
        ----------
        lon : ndarray,
            longitude of points to find
        lat : ndarray
            latitude of points to find
        grid : string, optional,
            "rho", "u", or "v" grid to search

        Returns
        -------
        indices : tuple of ndarray
            The indices for each dimension of the grid that are closest
            to the lon/lat points specified
        """

        glat = getattr(self, "lat_" + grid)
        glon = getattr(self, "lon_" + grid)
        xy = np.dstack([glat.ravel(), glon.ravel()])[0]
        pts = np.dstack([np.atleast_1d(lat), np.atleast_1d(lon)])[0]
        grid_tree = scipy.spatial.cKDTree(xy)
        dist, idx = grid_tree.query(pts)
        return np.unravel_index(idx, glat.shape)

    def ij(self, points):
        """
        Compute the fractional i,j indices of the grid from a
        set of lon, lat points.

        Parameters
        ----------
        points : list of tuples
            longitude, latitude points to compute i, j indicies

        Returns
        -------
        out : tuple of numpy masked array (with netcdf-type indexing),
            list of j,i indices for the given lon, lat points. NOTE: values
            that lie on the mask_rho are masked; however, if you wish to
            ignore masking, you can use the data field (i.data) directly.
            Values that do not lie within the grid are masked and stored as
            np.nan.

        Examples
        --------
        >>> a = ([-158, -160.5, -155.5], [20, 22.443, 19.5])
        >>> idx = g.ij(a)
        """

        from seapy.external.hindices import hindices

        # Interpolate the lat/lons onto the I, J
        xgrid, ygrid = np.ma.masked_equal(hindices(self.angle.T,
                                                   self.lon_rho.T, self.lat_rho.T,
                                                   points[0], points[1]), -999.0)
        mask = self.mask_rho[(ygrid.filled(0).astype(int),
                              xgrid.filled(0).astype(int))]
        xgrid[mask == 0] = np.ma.masked
        ygrid[mask == 0] = np.ma.masked
        return (ygrid, xgrid)

    def ijk(self, points, depth_adjust=False):
        """
        Compute the fractional i, j, k indices of the grid from a
        set of lon, lat, depth points.

        Parameters
        ----------
        points : list of tuples,
            longitude, latitude, depth points to compute i, j, k indicies.
            NOTE: depth is in meters (defaults to negative)
        depth_adjust : bool,
            If True, depths that are deeper (shallower) than the grid are set
            to the bottom (top) layer, 0 (N). If False, a nan value is used for
            values beyond the grid depth. Default is False.

        Returns
        -------
        out : tuple of numpy.maskedarray (with netcdf-type indexing),
            list of k, j, i indices for the given lon, lat, depth points

        Examples
        --------
        >>> a = ([-158, -160.5, -155.5], [20, 22.443, 19.5], [-10 -200 0])
        >>> idx = g.ijk(a)

        """
        # NOTE: Attempted to use a 3D griddata, but it took over 2 minutes
        # for each call, resulting in a 6minute runtime for this method
        # Reverted to 2D i,j indices, then looping a 1-D interpolation
        # to get depths for increased-speed (though this method is still slow)

        from scipy.interpolate import interp1d

        # Get the i,j points
        (j, i) = self.ij((points[0], points[1]))
        k = j * np.ma.masked
        grid_k = np.arange(0, self.n)
        depth = np.asanyarray(points[2])
        depth[depth > 0] *= -1

        # Determine the unique points
        good = np.where(~np.logical_or(i.mask, j.mask))[0]
        ii = np.floor(i[good]).astype(int)
        jj = np.floor(j[good]).astype(int)
        idx = seapy.unique_rows((jj, ii))
        fill_value = 0 if depth_adjust else np.nan
        for n in idx:
            pts = np.where(np.logical_and(jj == jj[n], ii == ii[n]))
            griddep = self.depth_rho[:, jj[n], ii[n]]
            if griddep[0] < griddep[-1]:
                griddep[-1] = 0.0
            else:
                griddep[0] = 0.0

            fi = interp1d(griddep, grid_k, bounds_error=False,
                          fill_value=fill_value)
            k[good[pts]] = fi(depth[good][pts])

        # Mask bad points
        l = np.isnan(k.data)
        i[l] = np.ma.masked
        j[l] = np.ma.masked
        k[l] = np.ma.masked

        return (k, j, i)

    def latlon(self, indices):
        """
        Compute the latitude and longitude from the given (i,j) indices
        of the grid

        Parameters
        ----------
        indices : list of tuples
            i, j points to compute latitude and longitude

        Returns
        -------
        out : tuple of ndarray
            list of lat,lon points from the given i,j indices

        Examples
        --------
        >>> a = [(23.4, 16.5), (3.66, 22.43)]
        >>> idx = g.latlon(a)
        """
        from scipy.interpolate import RegularGridInterpolator

        lati = RegularGridInterpolator((self.I[0, :], self.J[:, 0]),
                                       self.lat_rho.T)
        loni = RegularGridInterpolator((self.I[0, :], self.J[:, 0]),
                                       self.lon_rho.T)

        return (lati(indices), loni(indices))

    def rfactor(self):
        """
        Return the 2D field of the r-factor for the given grid.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray:
          array of r-factors size of the grid

        """
        hx = np.zeros(self.shape)
        hy = hx.copy()
        r = hx.copy()

        hx[:, :-1] = np.abs(np.diff(self.h, axis=1) /
                            (self.h[:, 1:] + self.h[:, :-1]))
        hy[:-1, :] = np.abs(np.diff(self.h, axis=0) /
                            (self.h[1:, :] + self.h[:-1, :]))
        hx[:, :-1] *= self.mask_u
        hy[:-1, :] *= self.mask_v

        r[:-1, :-1] = np.maximum(np.maximum(hx[:-1, :-1], hx[:-1, 1:]),
                                 np.maximum(hy[:-1, :-1], hy[1:, :-1]))
        r[:, -1] = r[:, -2]
        r[-1, :] = r[-2, :]
        hx = hy = 0
        return r * self.mask_rho

    def dHdxy(self):
        """
        Calculate the spatial derivative of water depth in each direction
        (xi and eta).

        Parameters
        ----------
        None

        Returns
        -------
        dHdxi : ndarray,
          Slope in x-direction
        dHdeta : ndarray,
          Slope in eta-direction
        """
        dHdxi = np.zeros(self.h.shape)
        dHdeta = np.zeros(self.h.shape)
        dHdxi[:, :-1] = -np.diff(self.h, axis=1) * self.pm[:, 1:]
        dHdxi[:, -1] = dHdxi[:, -2]
        dHdeta[:-1, :] = -np.diff(self.h, axis=0) * self.pn[1:, :]
        dHdeta[-1, :] = dHdeta[-2, :]

        return dHdxi, dHdeta

    def mask_poly(self, vertices, lat_lon=False, radius=0.0):
        """
        Create an np.masked_array of the same shape as the grid with values
        masked if they are not within the given polygon vertices

        Parameters
        ----------
        vertices: list of tuples,
            points that define the vertices of the polygon
        lat_lon : bool, optional,
            If True, the vertices are a list of lon, lat points rather
            than indexes

        Returns
        -------
        mask : np.masked_array
            mask of values that are located within the polygon

        Examples
        --------
        >>> vertices = [ (1,2), (4,5), (1,3) ]
        >>> mask = grid.mask_poly(vertices)
        """
        # If lat/lon vertices are given, we need to put these onto
        # the grid coordinates
        if lat_lon:
            points = self.ij(vertices, asint=True)
            vertices = list(zip(points[0], points[1]))

        # Now, with grid coordinates, test the grid against the vertices
        poly = matplotlib.path.Path(vertices)
        inside = poly.contains_points(np.vstack((self.J.flatten(),
                                                 self.I.flatten())).T,
                                      radius=radius)
        return np.ma.masked_where(inside.reshape(self.lat_rho.shape),
                                  np.ones(self.lat_rho.shape))
