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
    >>> m=seapy.mapping.map(region=(lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1]))
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
        region: tuple of latitude and longitude bounds
            (left, bottom, right, top)
        proj: string, optional
            projection to use (name from cartopy projection list, e.g.:
                 PlateCarree, LambertConformal, Miller, etc.)
        resolution: string, optional
            resolution to use for coastline, etc. From cartopy:
            'auto' (default), 'coarse', 'low', 'intermediate', 'high' or 'full'
        figsize: list, optional
            dimensions to use for creation of figure
        fig: matplotlib.pyplot.figure object, optional
            If you want to plot on a pre-configured figure, pass the figure object
            along with the axis object.
        fill_color: string, optional
            The color to use for the axis background
    """

    def __init__(self, grid=None, region=(-180, -40, 180, 40),
                 proj='PlateCarree', resolution='auto', figsize=(8., 6.),
                 fig=None, fill_color=cft.COLORS['water']):
        import inspect

        if grid is not None:
            grid = asgrid(grid)
            self.region = (np.min(grid.lon_rho), np.min(grid.lat_rho),
                           np.max(grid.lon_rho), np.max(grid.lat_rho))
        else:
            self.region = region

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
        self.clon = 0.5 * (region[0] + region[2])
        self.clat = 0.5 * (region[1] + region[3])

        args = inspect.getfullargspec(cproj)[0]
        if "central_latitude" not in args:
            self.clat = None

        self.fill = fill_color
        self.figsize = figsize
        self.res = resolution
        self.fig = fig
        if fig is None:
            self.new_figure()

    def new_figure(self, **kwargs):
        """
        Create a new mapping figure

        Parameters
        ----------
        fill_color: string, optional
           Color to fill the background of the axes with
        """
        self.pc = self.cb = self.ax = None

        # Create a figure
        self.fig = plt.figure(figsize=self.figsize, **kwargs)
        self.new_axes()

    def new_axes(self):
        """
        Create a new axes for plotting geospatial data
        """
        # Create the mapping axes
        if self.clat:
            self.ax = self.fig.add_subplot(1, 1, 1,
                                           projection=self.proj(
                                               central_longitude=self.clon,
                                               central_latitude=self.clat))
        else:
            self.ax = self.fig.add_subplot(1, 1, 1,
                                           projection=self.proj(
                                               central_longitude=self.clon))
        self.ax.set_extent((self.region[0], self.region[2],
                            self.region[1], self.region[3]),
                           crs=ccrs.PlateCarree())
        self.ax.set_facecolor(self.fill)

    def clear(self):
        """
        Clear the existing axes to draw a new frame
        """
        if self.ax:
            self.ax.clear()

    def land(self, facecolor="white", **kwargs):
        """
        Draw the land mask

        Parameters
        ----------
        facecolor: string, optional
            color to draw the mask with
        **kwargs: additional arguments to add_feature
        """
        self.ax.add_feature(cft.GSHHSFeature(self.res, [1]),
                            facecolor=facecolor, **kwargs)

    def gridlines(self, labels=(True, False, True, False), linewidth=1,
                  color='black', alpha=0.25, linestyle=":", **kwargs):
        """
        Draw grid lines

        Parameters
        ----------
        labels: boolean array
           True/False to display the grid labels on each side of the
           figure: [left, right, top, bottom]
        linewidth: size of grid lines
        color: color of grid lines
        alpha: transparency of grid lines
        linestyle: type of lines
        kwargs: additional arguments passed to the gridline
        """
        gl = self.ax.gridlines(draw_labels=True, linewidth=linewidth,
                               color=color, alpha=alpha,
                               linestyle=linestyle, **kwargs)
        for i, side in enumerate(("left", "right", "top", "bottom")):
            setattr(gl, f"{side}_labels", labels[i])

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
        self.pc = self.ax.pcolormesh(lon - dlon * 0.5, lat - dlat * 0.5,
                                     data, transform=self.proj(), **kwargs)

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
        self.pc = self.ax.contourf(lon, lat, data, transform=self.proj(),
                                   **kwargs)

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
        self.pc = self.ax.scatter(lon, lat, c=data, transform=self.proj(),
                                  **kwargs)

    def colorbar(self, label=None, location="bottom", **kwargs):
        """
        Display a colorbar on the figure

        Parameters
        ----------
        label: string, optional
            Colorbar label title
        location: string
            "bottom" (default), "top", "right", or "left"
        **kwargs: arguments, optional
            additional arguments to pass to colorbar
        """
        if self.pc is None:
            return

        l, b, w, h = self.ax.get_position().bounds

        location = location.lower()
        if location == "right":
            cax = self.fig.add_axes([0.94, 0.25, 0.03, h * 0.8])
            orientation = "vertical"
        elif location == "left":
            cax = self.fig.add_axes([0.01, 0.25, 0.03, h * 0.8])
            orientation = "vertical"
        elif location == "bottom":
            cax = self.fig.add_axes([0.25, b - 0.07, 0.5, 0.03])
            orientation = "horizontal"
        elif location == "top":
            cax = self.fig.add_axes([0.25, b + 1.13 * h, 0.5, 0.03])
            orientation = "horizontal"

        self.cb = plt.colorbar(self.pc, cax=cax, orientation=orientation,
                               shrink=0.8, **kwargs)
        # self.cb = plt.colorbar(self.pc, ax=self.ax, location=location,
        #                        **kwargs)
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

    def __init__(self, grid=None, region=(-163, 17, -153, 24),
                 figsize=(8., 6.)):
        super().__init__(grid=grid, region=region, figsize=figsize)

    def land(self, facecolor="white", **kwargs):
        """
        Draw the GIS coastline data from the state of Hawaii to draw the
        land boundaries. This does not include rivers, etc., only the
        coastline.

        Parameters
        ----------
        color: string, optional
            Color to draw the land mask with
        **kwargs:
            Additional arguments to add_feature

        Returns
        -------
        None

        """
        _shape_file = os.path.dirname(__file__) + "/hawaii_coast/hawaii.shp"

        if not hasattr(self, "feature"):
            self.feature = cft.ShapelyFeature(
                cartopy.io.shapereader.Reader(_shape_file).geometries(),
                self.proj())

        # Draw the loaded shapes
        self.ax.add_feature(self.feature, facecolor=facecolor, **kwargs)
