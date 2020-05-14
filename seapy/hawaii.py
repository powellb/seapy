#!/usr/bin/env python
"""
  hawaii.py

  State Estimation and Analysis for PYthon

  Utilities for dealing with data around Hawaii

    Examples
    --------

    Assume you have longitude, latitude, and sst values:

    >>> m=seapy.hawaii()
    >>> m.pcolormesh(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
    >>> m.land()
    >>> m.colorbar(label="Sea Surface Temp [$^\circ$C]",cticks=[22,23,24,25,26])
    >>> m.ax.patch.set_facecolor("aqua")
    >>> m.ax.patch.set_alpha(1)
    >>> m.fig.patch.set_alpha(0.0)
    >>> m.fig.savefig("sst.png",dpi=100)

  Written by Brian Powell on 9/4/14
  Copyright (c)2020 University of Hawaii under the MIT-License.
"""


from .mapping import map
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import os

_shape_file = os.path.dirname(__file__) + "/hawaii_coast/hawaii"


class hawaii(map):

    def __init__(self, grid=None, llcrnrlon=-163, llcrnrlat=17, urcrnrlon=-153,
                 urcrnrlat=24, figsize=(8., 6.), dlat=1, dlon=2, fig=None, ax=None,
                 fill_color="aqua"):
        super().__init__(grid=grid, llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                         urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                         figsize=figsize, dlat=dlat, dlon=dlon, fig=fig, ax=ax,
                         fill_color=fill_color)

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

        if hasattr(self.basemap, "coast") == False or hasattr(self, "landpoly"):
            self.basemap.readshapefile(_shape_file, "coast")
            vert = []
            for shape in self.basemap.coast:
                vert.append(shape)

            self.landpoly = PolyCollection(
                vert, facecolors=color, edgecolors=color)
        # Draw the loaded shapes
        self.ax.add_collection(self.landpoly)
