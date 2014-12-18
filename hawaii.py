#!/usr/bin/env python
"""
  hawaii.py
  
  State Estimation and Analysis for PYthon

  Utilities for dealing with data around Hawaii

  Written by Brian Powell on 9/4/14
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

from .map import map
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import os

_shape_file = os.path.dirname(__file__)+"/hawaii_coast/hawaii"

class hawaii(map):
    def __init__(self, figsize=(8.,6.), dlat=1, dlon=2):
        super().__init__(llcrnrlon=-164.5, llcrnrlat=16.6,
                   urcrnrlon=-152.0, urcrnrlat=24.2, figsize=figsize, 
                   dlat=dlat, dlon=dlon)

    def land(self, color="black"):
        """
        Use the GIS coastline data from the state of Hawaii to draw the
        land boundaries. This does not delinite rivers, etc., only the
        coastline.
        """
        
        if  hasattr(self.basemap,"coast") == False or hasattr(self, "landpoly"):
            self.basemap.readshapefile(_shape_file, "coast")
            vert=[]
            for shape in self.basemap.coast:
                vert.append(shape)

            self.landpoly = PolyCollection(vert,facecolors=color,edgecolors=color)
        # Draw the loaded shapes
        self.ax.add_collection(self.landpoly)
