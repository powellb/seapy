#!/usr/bin/env python
"""
  hawaii.py
  
  State Estimation and Analysis for PYthon

  Utilities for dealing with data around Hawaii

  Written by Brian Powell on 9/4/14
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
from matplotlib.collections import PolyCollection
import os

_shape_file = os.path.dirname(__file__)+"/hawaii_coast/hawaii"

class map(object):
    def __init__(self, figsize=(8.,6.), color="black", dlat=1, dlon=2):
        self.basemap = Basemap(llcrnrlon=-164.5, llcrnrlat=16.7,
                               urcrnrlon=-152.25, urcrnrlat=24,
                               projection='lcc', lat_0=20,
                               lon_0=-158.0, resolution='c', area_thresh=0.0)
                         
        self.fig = plt.figure(figsize=figsize)
        self.ax = self.fig.add_axes([-0.01, 0.27, 1.01, 0.7])
        self.basemap.drawmapboundary(fill_color="aqua")

        self.color=color

        # Create the lat/lon lines
        self.basemap.drawmeridians(np.arange(-164.5,-152.25,dlon),color="0.5",
            linewidth=0.25, dashes=[1,1,0,1,0], labels=[0,0,0,1],fontsize=12)
        self.basemap.drawparallels(np.arange(16.7,24,dlat),color="0.5",
            linewidth=0.25, dashes=[1,1,0,1,0], labels=[1,0,0,0],fontsize=12)
            
    def pcolor(self, lon, lat, data, label=None, cticks=None, **kwargs):
        # Pcolor requires a modification to the locations to line up with
        # the geography
        dlon=lon*0;
        dlat=lat*0;
        dlon[:,0:-1]=lon[:,1:]-lon[:,0:-1]
        dlat[0:-1,:]=lat[1:,:]-lat[0:-1,:]
        x,y = self.basemap(lon-dlon*0.5,lat-dlat*0.5)

        self.pc = self.basemap.pcolor(x,y,data,**kwargs)
        self.cax = self.fig.add_axes([0.25, 0.18, 0.5, 0.03])
        self.cb = plt.colorbar(self.pc, cax=self.cax, orientation="horizontal",
                                ticks=cticks)
        self.basemap.set_axes_limits(ax=self.ax)
        if label != None:
            self.cb.set_label(label)
        self.land()
    
    def land(self):
        if self.color != None:
            self.basemap.readshapefile(_shape_file, "coast", drawbounds=False)
            vert=[]
            for shape in self.basemap.coast:
                vert.append(shape)

            self.land = PolyCollection(vert,facecolors=self.color, 
                               edgecolors=self.color)
            self.ax.add_collection(self.land)
        else:
            self.basemap.readshapefile(_shape_file, "coast", drawbounds=True)

        