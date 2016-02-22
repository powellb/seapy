#!/usr/bin/env python
"""
  map.py

  State Estimation and Analysis for PYthon

  Utilities for dealing with basemap plotting. These routnes are simply
  abstractions over the existing basemap to make it quicker for generating
  basemap plots and figures.

    Examples
    -------

    Assume you have longitude, latitude, and sst values:

    >>> m=seapy.mapping.map(llcrnrlon=lon[0,0],llcrnrlat=lat[0,0],
    >>>     urcrnrlon=lon[-1,-1],urcrnrlat=lat[-1,-1],dlat=2,dlon=2)
    >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
    >>> m.land()
    >>> m.colorbar(label="Sea Surface Temp [$^\circ$C]",cticks=[22,23,24,25,26])
    >>> m.ax.patch.set_facecolor("aqua")
    >>> m.ax.patch.set_alpha(1)
    >>> m.fig.patch.set_alpha(0.0)
    >>> m.fig.savefig("sst.png",dpi=100)


  Written by Brian Powell on 9/4/14
  Copyright (c)2016 University of Hawaii under the BSD-License.
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
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

    def __init__(self, grid=None, llcrnrlon=-180, llcrnrlat=-40, urcrnrlon=180,
                 urcrnrlat=40, proj='lcc', resolution='c', figsize=(8., 6.),
                 dlat=1, dlon=2):
        """
        map class for abstracting the basemap methods for quick and easy creation
        of geographically referenced data figures


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
            projection to use for map
        resolution: character
            resolution to use for coastline, etc. From Basemap:
            'c' (crude), 'l' (low), 'i' (intermediate),
            'h' (high), 'f' (full), or None
        figsize: list, optional
            dimensions to use for creation of figure
        dlat: float, optional
            how often to mark latitude lines
        dlon: float, optional
            how often to mark longitude lines


        Returns
        -------
        None

        """
        if grid is not None:
            grid = asgrid(grid)
            llcrnrlat = np.min(grid.lat_rho)
            urcrnrlat = np.max(grid.lat_rho)
            llcrnrlon = np.min(grid.lon_rho)
            urcrnrlon = np.max(grid.lon_rho)

        self.basemap = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                               urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                               projection=proj,
                               lat_0=urcrnrlat - (urcrnrlat - llcrnrlat) / 2.,
                               lon_0=urcrnrlon - (urcrnrlon - llcrnrlon) / 2.,
                               resolution=resolution, area_thresh=0.0)

        self.figsize = figsize
        self.dlon = dlon
        self.dlat = dlat
        self.fig = None
        self.new_figure()

    def new_figure(self, fill_color="aqua"):
        """
        Create a new figure for plotting
        """
        if self.fig is not None:
            self.ax.set_axis_off()
            plt.close(self.fig)

        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_axes([-0.01, 0.25, 1.01, 0.7])
        self.basemap.drawmapboundary(fill_color=fill_color)
        # Create the lat/lon lines
        delta = self.basemap.urcrnrlon - self.basemap.llcrnrlon
        nticks = int(delta / self.dlon)
        if delta / nticks > 1:
            lon_lines = np.linspace(int(self.basemap.llcrnrlon - self.dlon),
                                    int(self.basemap.urcrnrlon + self.dlon), nticks + 2)
        else:
            lon_lines = np.linspace(self.basemap.llcrnrlon - self.dlon,
                                    self.basemap.urcrnrlon + self.dlon, nticks + 2)
        # lon_lines = np.arange(self.basemap.llcrnrlon,
        #     self.basemap.urcrnrlon, self.dlon)
        self.basemap.drawmeridians(lon_lines, color="0.5",
                                   linewidth=0.25, dashes=[1, 1, 0.1, 1], labels=[0, 0, 0, 1], fontsize=12)
        delta = self.basemap.urcrnrlat - self.basemap.llcrnrlat
        nticks = int(delta / self.dlat)
        if delta / nticks > 1:
            lat_lines = np.linspace(int(self.basemap.llcrnrlat - self.dlat),
                                    int(self.basemap.urcrnrlat + self.dlat), nticks + 2)
        else:
            lat_lines = np.linspace(self.basemap.llcrnrlat - self.dlat,
                                    self.basemap.urcrnrlat + self.dlat, nticks + 2)
        # lat_lines = np.arange(self.basemap.llcrnrlat,
        #     self.basemap.urcrnrlat, self.dlat)
        self.basemap.drawparallels(lat_lines, color="0.5",
                                   linewidth=0.25, dashes=[1, 1, 0.1, 1], labels=[1, 0, 0, 0], fontsize=12)

    def land(self, color="black"):
        """
        Draw the land mask

        Parameters
        ----------
        color: string, optional
            color to draw the mask with
        """
        self.basemap.drawcoastlines()
        self.basemap.drawcountries()
        self.basemap.fillcontinents(color=color)

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
        self.ax.set_xlim(x)
        self.ax.set_ylim(y)
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
        x, y = self.basemap(lon - dlon * 0.5, lat - dlat * 0.5)
        self.pc = self.ax.pcolormesh(x, y, data, **kwargs)

    def colorbar(self, label=None, cticks=None, **kwargs):
        """
        Display a colorbar on the figure

        Parameters
        ----------
        label: string, optional
            Colorbar label title
        cticks: array, optional
            Where to place the tick marks and values for the colorbar
        """
        self.cax = self.fig.add_axes([0.25, 0.16, 0.5, 0.03])
        self.cb = plt.colorbar(self.pc, cax=self.cax, orientation="horizontal",
                               ticks=cticks)
        self.basemap.set_axes_limits(ax=self.ax)
        if label is not None:
            self.cb.set_label(label)


