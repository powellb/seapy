#!/usr/bin/env python
"""
  mapping.py

  State Estimation and Analysis for PYthon

  Utilities for dealing with plotting maps using cartopy. These routines are
  provide simplified abstractions over the existing cartopy to make it quicker
  for generating plots and figures.

  Copyright (c)2010--2025 University of Hawaii under the MIT-License.

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
    Create a new mapping object for producing geographic figures.

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

    Examples
    --------
    >>> m=seapy.mapping.map(region=(lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1]))
    >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
    >>> m.land()
    >>> m.colorbar(label="Sea Surface Temp",cticks=[22,23,24,25,26])
    >>> m.ax.patch.set_alpha(1)
    >>> m.fig.patch.set_alpha(0.0)
    >>> m.fig.savefig("sst.png",dpi=100)

"""

    def __init__(self, nrows=1, ncols=1, grid=None, region=(-180, -40, 180, 40),
                 proj='PlateCarree', resolution='auto', figsize=(8., 6.),
                 fig=None, fill_color=cft.COLORS['water'], **kwargs):

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

        self.nrows = nrows
        self.ncols = ncols
        self.fill = fill_color
        self.figsize = figsize
        self.res = resolution
        self.fig = fig
        if fig is None:
            self.new_figure(**kwargs)

    def new_figure(self, **kwargs):
        """
        Create a new mapping figure

        Parameters
        ----------
        fill_color: string, optional
           Color to fill the background of the axes with
        """
        # Create the mapping setup
        self.cur_sp = self.cur_pc = None

        if self.clat:
            geokw = dict(projection=self.proj(central_longitude=self.clon,
                                              central_latitude=self.clat))
        else:
            geokw = dict(projection=self.proj(central_longitude=self.clon))

        # Create a figure
        self.fig, self.ax = plt.subplots(self.nrows, self.ncols,
                                         figsize=self.figsize, **kwargs,
                                         subplot_kw=geokw)
        self.ax = np.atleast_2d(self.ax)
        self.pc = np.empty(self.ax.shape, dtype=object)

        for ax in self.fig.axes:
            ax.set_extent((self.region[0], self.region[2],
                           self.region[1], self.region[3]),
                          crs=ccrs.PlateCarree())
            ax.set_facecolor(self.fill)

        # Set the Current Axis
        self.subplot()

    def subplot(self, ax=None):
        """
        Set the current subplot axis to work on. Subplots start in the upper,
        left with (0, 0) down to the lower right with (nrows-1, ncolumns-1).

        Parameters
        ----------
        ax : tuple
          A tuple of the (column, row) axis that you wish to currently
          work with

        Examples
        --------
        >>> m=seapy.mapping.map(region=(lon[0,0],lat[0,0],lon[-1,-1],lat[-1,-1]))
        >>> m.plot(lon, lat, 'k')
        >>> m.subplot((0,1))
        >>> m.pcolormesh(lon, lat, data)

        """
        self.cur_ax = ax if ax else (0, 0)
        self.cur_sp = self.ax[self.cur_ax]

    def clear(self):
        """
        Clear the existing axes to draw a new frame
        """
        self.cur_sp.clear()

    def title(self, *args, **kwargs):
        """
        Set the title of the current plot

        Parameters
        ----------
        Same as matplotlib.axes.Axes.set_title
        """
        self.cur_sp.set_title(*args, **kwargs)

    def land(self, facecolor="white", **kwargs):
        """
        Draw the land mask

        Parameters
        ----------
        facecolor: string, optional
            color to draw the mask with. If you specify 'none', it will not
            fill land.
        **kwargs: additional arguments to add_feature
        """
        self.cur_sp.add_feature(cft.GSHHSFeature(self.res, [1]),
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
        gl = self.cur_sp.gridlines(draw_labels=True, linewidth=linewidth,
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
            additional arguments to pass to pcolormesh
        """
        # Pcolor requires a modification to the locations to line up with
        # the geography
        dlon = lon * 0
        dlat = lat * 0
        dlon[:, 0:-1] = lon[:, 1:] - lon[:, 0:-1]
        dlat[0:-1, :] = lat[1:, :] - lat[0:-1, :]
        self.pc[self.cur_ax] = self.cur_sp.pcolormesh(lon - dlon * 0.5,
                                                      lat - dlat * 0.5,
                                                      data, transform=self.proj(),
                                                      **kwargs)

    def text(self, lon, lat, text, **kwargs):
        """
        plot data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        text: string
            Text to display
        **kwargs: arguments, optional
            additional arguments to pass to plot
        """
        self.cur_sp.text(lon, lat, text, transform=self.proj(), **kwargs)

    def plot(self, lon, lat, *args, **kwargs):
        """
        plot data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        **kwargs: arguments, optional
            additional arguments to pass to plot
        """
        return self.cur_sp.plot(lon, lat, *args, transform=self.proj(), **kwargs)

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
            additional arguments to pass to contourf
        """
        self.pc[self.cur_ax] = self.cur_sp.contourf(lon, lat, data,
                                                    transform=self.proj(),
                                                    **kwargs)

    def contour(self, lon, lat, data, **kwargs):
        """
        contour field data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        data: array
            data to contourf
        **kwargs: arguments, optional
            additional arguments to pass to contourf
        """
        self.pc[self.cur_ax] = self.cur_sp.contour(lon, lat, data,
                                                   transform=self.proj(),
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
            data to scatter
        **kwargs: arguments, optional
            additional arguments to pass to scatter
        """
        self.pc[self.cur_ax] = self.cur_sp.scatter(lon, lat, c=data,
                                                   transform=self.proj(),
                                                   **kwargs)

    def streamplot(self, lon, lat, u, v, **kwargs):
        """
        streamline plot vector data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        u: array
            u-component to display
        v: array
            v-component to display
        **kwargs: arguments, optional
            additional arguments to pass to streamplot
        """
        return self.cur_sp.streamplot(lon, lat, u, v, transform=self.proj(),
                                      **kwargs)

    def quiver(self, lon, lat, u, v, **kwargs):
        """
        quiver plot vector data onto our geographic plot

        Parameters
        ----------
        lon: array
            Longitude field for data
        lat: array
            Latitude field for data
        u: array
            u-component to display
        v: array
            v-component to display
        **kwargs: arguments, optional
            additional arguments to pass to quiver
        """
        return self.cur_sp.quiver(lon, lat, u, v, transform=self.proj(),
                                  **kwargs)

    def colorbar(self, ax=None, label=None, location="bottom", newaxis=False,
                 pad=0.1, shrink=0.8, **kwargs):
        """
        Display a colorbar. The colorbar can be attached to
        individual subplots [default] or it can span multiple
        subplots. If it spans multiple, it will use the values
        from the first axis of the span for the colorbar range.

        Parameters
        ----------
        ax: index or slice, optional
            The axis (or more) that the colorbar should span. The default
            is the current axis.
        label: string, optional
            Colorbar label title
        location: string
            "bottom" (default), "top", "right", or "left"
        newaxis : bool
            Create a new axis outside of the axes. Results will vary
            depending on the way the figure is constructed, so you
            must know what you are doing. This is only for creating
            colorbars outside of the subplots.
        pad : float
            How large of a pad to space the colorbar if creating a new
            axis.
        shrink : float
            The percentage to shrink the colorbar by
        **kwargs: arguments, optional
            additional arguments to pass to colorbar

        Returns
        -------
        cb : matplotlib.colorbar.Colorbar
            the colorbar object that was created (only needed if you
            wish to alter the colorbar)
        """

        # Set the values object(s) and axes object(s)
        pc = self.pc[ax].flatten()[0] if ax else self.pc[self.cur_ax]
        ax = self.ax[ax] if ax else self.cur_sp

        # If user wants a new axis, calculate a new axis to hold
        # the colorbar.
        if newaxis or len(self.ax) == 1:
            ax = np.atleast_1d(ax)
            sz = np.zeros((ax.size, 4))
            for i, a in enumerate(ax.flatten()):
                sz[i, :] = np.array(a.get_position().bounds)

            if location == "right":
                orientation = "vertical"
                ll = np.max(sz[:, 0] + sz[:, 2]) + pad
                _, idx = np.unique(sz[:, 1], return_index=True)
                hh = np.minimum(0.8, np.sum(sz[idx, 3]))
                ww = 0.03
                bb = np.min(sz[:, 1]) + hh / 2 - hh * shrink / 2
                hh *= shrink

            elif location == "left":
                orientation = "vertical"
                ll = np.min(sz[:, 0]) - pad * 2
                _, idx = np.unique(sz[:, 1], return_index=True)
                hh = np.minimum(0.8, np.sum(sz[idx, 3]))
                ww = 0.03
                bb = np.min(sz[:, 1]) + hh / 2 - hh * shrink / 2
                hh *= shrink

            elif location == "top":
                orientation = "horizontal"
                bb = np.max(sz[:, 1] + sz[:, 3]) + pad * 2
                _, idx = np.unique(sz[:, 0], return_index=True)
                ww = np.minimum(0.8, np.sum(sz[idx, 2]))
                hh = 0.03
                ll = np.min(sz[:, 0]) + ww / 2 - ww * shrink / 2
                ww *= shrink

            else:  # bottom
                orientation = "horizontal"
                bb = np.min(sz[:, 1]) - pad * 2
                _, idx = np.unique(sz[:, 0], return_index=True)
                ww = np.minimum(0.8, np.sum(sz[idx, 2]))
                hh = 0.03
                ll = np.min(sz[:, 0]) + ww / 2 - ww * shrink / 2
                ww *= shrink

            cax = self.fig.add_axes([ll, bb, ww, hh])
            cb = self.fig.colorbar(pc, cax=cax, orientation=orientation,
                                   shrink=shrink)
        else:
            # Let matplotlib do its own colorbar adjustment
            cb = self.fig.colorbar(pc, ax=ax, location=location,
                                   shrink=shrink, **kwargs)
        if label is not None:
            cb.set_label(label)

        return cb


class hawaii(map):
    """
      Make plots using high resolution coastlines around Hawaii

        Examples
        --------
        >>> m=seapy.hawaii()
        >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
        >>> m.land()
        >>> m.colorbar(label="Sea Surface Temp",cticks=[22,23,24,25,26])
        >>> m.ax.patch.set_alpha(1)
        >>> m.fig.patch.set_alpha(0.0)
        >>> m.fig.savefig("sst.png",dpi=100)

      Copyright (c)2010--2025 University of Hawaii under the MIT-License.
    """
    # Class variable for the shapes
    import importlib.resources as importer
    _shape_file = importer.files("seapy").joinpath("hawaii_coast", "hawaii.shp")

    def __init__(self, *args, region=(-163, 17, -153, 24), figsize=(8., 6.),
                 **kwargs):
        super(hawaii, self).__init__(*args, region=region, figsize=figsize, **kwargs)

        self.shape = cft.ShapelyFeature(
            cartopy.io.shapereader.Reader(self._shape_file).geometries(),
            self.proj())

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

        # Draw the loaded shapes
        self.cur_sp.add_feature(self.shape, facecolor=facecolor, **kwargs)
