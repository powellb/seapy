#!/usr/bin/env python
"""
  mapping.py

  State Estimation and Analysis for PYthon

  Utilities for dealing with plotting maps using cartopy. These routines are
  provide simplified abstractions over the existing cartopy to make it quicker
  for generating plots and figures.

  Copyright (c)2010--2021 University of Hawaii under the MIT-License.

"""
import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cft
import os
from warnings import warn
from seapy.model import asgrid


def gen_coastline(lon, lat, bathy, depth=0):
    """
    Given lon, lat, and bathymetry, generate vectors of line segments
    of the coastline. This can be exported to matlab (via savemat) to be
    used with the 'editmask' routine for creating grid masks.

    Input
    -----
    lon : array,
        longitudes of bathymetry locations
    lat : array,
        latitudes of bathymetry locations
    bathy : array,
        bathymetry (negative for ocean, positive for land) values
    depth : float,
        depth to use as the definition of the coast

    Returns
    -------
    lon : ndarray,
        vector of coastlines, separated by nan (matlab-style)
    lat : ndarray,
        vector of coastlines, separated by nan (matlab-style)
    """
    CS = plt.contour(lon, lat, bathy, [depth - 0.25, depth + 0.25])
    lon = list()
    lat = list()
    for col in CS.collections:
        for path in col.get_paths():
            lon.append(path.vertices[:, 0])
            lon.append(np.nan)
            lat.append(path.vertices[:, 1])
            lat.append(np.nan)
    return (np.hstack(lon), np.hstack(lat))


class map(object):
    """
    Examples
    -------
    >>> m=seapy.mapping.map(llcrnrlon=lon[0,0],llcrnrlat=lat[0,0],
    >>>     urcrnrlon=lon[-1,-1],urcrnrlat=lat[-1,-1],dlat=2,dlon=2)
    >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
    >>> m.land()
    >>> m.colorbar(label="Sea Surface Temp [$^\circ$C]",cticks=[22,23,24,25,26])
    >>> m.ax.patch.set_alpha(1)
    >>> m.fig.patch.set_alpha(0.0)
    >>> m.fig.savefig("sst.png",dpi=100)

    Parameters
    ----------
        grid: seapy.model.grid or string, optional:
            grid to use to define boundaries
        llcrnrlon: float, optional
            longitude of lower, left corner
        llcrnrlat: float, optional
            latitude of lower, left corner
        urcrnrlon: float, optional
            longitude of upper, right corner
        urcrnrlat: float, optional
            latitude of upper, right corner
        proj: string, optional
            projection to use (name from cartopy projection list, e.g.:
                 PlateCarree, LambertConformal, Miller, etc.)
        resolution: string, optional
            resolution to use for coastline, etc. From cartopy:
            'auto' (default), 'coarse', 'low', 'intermediate', 'high' or 'full'
        figsize: list, optional
            dimensions to use for creation of figure
        dlat: float, optional
            interval to mark latitude lines (e.g., if dlat=0.5 every 0.5deg mark)
        dlon: float, optional
            interval to mark longitude lines (e.g., if dlon=0.5 every 0.5deg mark)
        fig: matplotlib.pyplot.figure object, optional
            If you want to plot on a pre-configured figure, pass the figure object
            along with the axis object.
        ax: matplotlib.pyplot.axis object, optional
            If you want to plot on a pre-configured figure, pass the axis object
            along with the figure object.
        fill_color: string, optional
            The color to use for the axis background
    """

    def __init__(self, grid=None, llcrnrlon=-180, llcrnrlat=-40,
                 urcrnrlon=180, urcrnrlat=40,
                 proj='PlateCarree',
                 resolution='auto', figsize=(8., 6.),
                 dlat=1, dlon=2, fig=None, ax=None,
                 fontsize=12, fill_color=cft.COLORS['water']):
        import inspect

        if grid is not None:
            grid = asgrid(grid)
            llcrnrlat = np.min(grid.lat_rho)
            urcrnrlat = np.max(grid.lat_rho)
            llcrnrlon = np.min(grid.lon_rho)
            urcrnrlon = np.max(grid.lon_rho)

        # Set up the projection scheme
        try:
            cproj = getattr(ccrs, proj)
        except AttributeError:
            warn(
                f"WARNING: {proj} is not a valid projection. " +
                "Using 'PlateCarree' instead")
            cproj = ccrs.PlateCarree
        self.proj = cproj

        # Determine what our projection needs
        clon = 0.5 * (urcrnrlon + llcrnrlon)
        clat = 0.5 * (urcrnrlat + llcrnrlat)
        self.fig = plt.figure(figsize=figsize)
        args = inspect.getfullargspec(cproj)[0]
        if "central_latitude" in args:
            self.ax = self.fig.add_subplot(1, 1, 1,
                                           projection=cproj(central_longitude=clon,
                                                            central_latitude=clat))
        else:
            self.ax = self.fig.add_subplot(1, 1, 1,
                                           projection=cproj(central_longitude=clon))
        self.ax.set_extent((llcrnrlon, urcrnrlon, llcrnrlat, urcrnrlat),
                           crs=ccrs.PlateCarree())
        self.ax.set_facecolor(fill_color)
        self.fill = fill_color
        self.figsize = figsize
        delta = (np.abs(urcrnrlon - llcrnrlon) // .04) / 100
        self.dlon = np.minimum(delta, dlon)
        delta = (np.abs(urcrnrlat - llcrnrlat) // .04) / 100
        self.dlat = np.minimum(delta, dlat)
        self.res = resolution
        self.fontsize = fontsize
        # reset = True if fig is None else False
        # self.new_figure(reset=reset)

    def new_figure(self, fill_color=None, reset=False, dpi=150):
        """
        Create or update a figure for plotting

        Parameters
        ----------
        fill_color: string, optional
           Color to fill the background of the axes with
        reset: bool, optional
           Reset the figure
        """
        if reset:
            if self.ax:
                self.ax.set_axis_off()
                self.ax = None
            if self.fig:
                self.fig.clf()
                self.fig = None

        if self.fig is None or self.ax is None:
            self.fig = plt.figure(figsize=self.figsize, dpi=dpi)
            self.ax = self.fig.add_axes([-0.01, 0.25, 1.01, 0.7])

        if fill_color is None:
            fill_color = self.fill_color

        # Create the longitude lines

        # Create the latitude lines

    def land(self, color="black"):
        """
        Draw the land mask

        Parameters
        ----------
        color: string, optional
            color to draw the mask with
        """
        self.ax.add_feature(cft.GSHHSFeature(self.res, [1]),
                            facecolor=color)

    def zoom(self, xrange, yrange):
        """
        zoom the figure to a specified lat, lon range

        Parameters
        ----------
        xrange: array
            minimum and maximum longitudes to display
        yrange: array
            minimum and maximum latitudes to display
        """
        x, y = self.basemap(xrange, yrange)
        self.ax.set_limits((xrange.min(), xrange.max(),
                            yrange.min(), yrange.max()))
        self.fig.canvas.draw()

    def pcolormesh(self, lon, lat, data, **kwargs):
        """
        pcolormesh field data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        data: array
            data to pcolor
        **kwargs: arguments, optional
            additional arguments to pass to pcolor
        """
        # Pcolor requires a modification to the locations to line up with
        # the geography
        dlon = lon * 0
        dlat = lat * 0
        dlon[:, 0:-1] = lon[:, 1:] - lon[:, 0:-1]
        dlat[0:-1, :] = lat[1:, :] - lat[0:-1, :]
        self.pc = self.ax.pcolormesh(lon - dlon * 0.5, lon - dlon * 0.5,
                                     data, **kwargs)

    def contourf(self, lon, lat, data, **kwargs):
        """
        contourf field data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        data: array
            data to contourf
        **kwargs: arguments, optional
            additional arguments to pass to pcolor
        """
        self.pc = self.ax.contourf(lon, lat, data, **kwargs)

    def scatter(self, lon, lat, data, **kwargs):
        """
        scatter plot data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        data: array
            data to pcolor
        **kwargs: arguments, optional
            additional arguments to pass to pcolor
        """
        self.pc = self.ax.scatter(lon, lat, c=data, **kwargs)

    def colorbar(self, label=None, cticks=None, **kwargs):
        """
        Display a colorbar on the figure

        Parameters
        ----------
        label: string, optional
            Colorbar label title
        cticks: array, optional
            Where to place the tick marks and values for the colorbar
        **kwargs: arguments, optional
            additional arguments to pass to colorbar
        """
        self.cax = self.fig.add_axes([0.25, 0.16, 0.5, 0.03])
        self.cb = plt.colorbar(self.pc, cax=self.cax, orientation="horizontal",
                               ticks=cticks, **kwargs)
        if label is not None:
            self.cb.set_label(label)


class hawaii(map):
    """
      Make plots using high resolution coastlines around Hawaii

        Examples
        --------
        >>> m=seapy.hawaii()
        >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
        >>> m.land()
        >>> m.colorbar(label="Sea Surface Temp [$^\circ$C]",cticks=[22,23,24,25,26])
        >>> m.ax.patch.set_alpha(1)
        >>> m.fig.patch.set_alpha(0.0)
        >>> m.fig.savefig("sst.png",dpi=100)

      Copyright (c)2010--2021 University of Hawaii under the MIT-License.
    """

    def __init__(self, grid=None, llcrnrlon=-163, llcrnrlat=17,
                 urcrnrlon=-153, urcrnrlat=24, figsize=(8., 6.),
                 dlat=1, dlon=2):
        super().__init__(grid=grid, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                         urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                         figsize=figsize, dlat=dlat, dlon=dlon)

    def land(self, color="black"):
        """
        Draw the GIS coastline data from the state of Hawaii to draw the
        land boundaries. This does not include rivers, etc., only the
        coastline.

        Parameters
        ----------
        color: string, optional
            Color to draw the land mask with

        Returns
        -------
        None

        """
        _shape_file = os.path.dirname(__file__) + "/hawaii_coast/hawaii.shp"

        if not hasattr(self, "feature"):
            print('load data')
            self.feature = cft.ShapelyFeature(
                cartopy.io.shapereader.Reader(_shape_file).geometries(),
                self.proj(central_longitude=-158))

        # Draw the loaded shapes
        self.ax.add_feature(self.feature, facecolor=color)
