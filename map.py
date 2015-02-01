#!/usr/bin/env python
"""
  map.py
  
  State Estimation and Analysis for PYthon

  Utilities for dealing with basemap plotting. These routnes are simply
  abstractions over the existing basemap to make it quicker for generating
  basemap plots and figures.

    Example
    -------
    
    Assume you have longitude, latitude, and sst values:
    
    m=seapy.map(llcrnrlon=lon[0,0],llcrnrlat=lat[0,0],
        urcrnrlon=lon[-1,-1],urcrnrlat=lat[-1,-1],dlat=2,dlon=2)
    m.pcolor(lon,lat,sst,vmin=22,vmax=26,cmap=plt.cm.bwr)
    m.land()
    m.colorbar(label="Sea Surface Temp [$^\circ$C]",cticks=[22,23,24,25,26])
    m.ax.patch.set_facecolor("aqua")
    m.ax.patch.set_alpha(1)
    m.fig.patch.set_alpha(0.0)
    m.fig.savefig("sst.png",dpi=100)
    
    
  Written by Brian Powell on 9/4/14
  Copyright (c)2013 University of Hawaii under the BSD-License.
"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os

class map(object):
    def __init__(self, llcrnrlon=-180, llcrnrlat=-40, urcrnrlon=180, 
                 urcrnrlat=40, figsize=(8.,6.), dlat=1, dlon=2):
        """
        map class for abstracting the basemap methods for quick and easy creation
        of geographically referenced data figures
    
    
        Parameters
        ----------
        llcrnrlon: float, optional
            longitude of lower, left corner
        llcrnrlat: float, optional
            latitude of lower, left corner
        urcrnrlon: float, optional
            longitude of upper, right corner
        urcrnrlat: float, optional
            latitude of upper, right corner
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
        self.basemap = Basemap(llcrnrlon=llcrnrlon, llcrnrlat=llcrnrlat,
                               urcrnrlon=urcrnrlon, urcrnrlat=urcrnrlat,
                               projection='lcc', 
                               lat_0=urcrnrlat-(urcrnrlat-llcrnrlat)/2.,
                               lon_0=urcrnrlon-(urcrnrlon-llcrnrlon)/2.,
                               resolution='c', area_thresh=0.0)
        self.figsize=figsize
        self.dlon=dlon
        self.dlat=dlat
        self.fig=None
        self.new_figure()

    def new_figure(self):
        """
        Create a new figure for plotting
        """
        if self.fig != None:
            self.ax.set_axis_off()
            plt.close(self.fig)
            
        self.fig = plt.figure(figsize=self.figsize)
        self.ax = self.fig.add_axes([-0.01, 0.25, 1.01, 0.7])
        self.basemap.drawmapboundary(fill_color="aqua")
        # Create the lat/lon lines
        self.basemap.drawmeridians(np.arange(self.basemap.llcrnrlon,
            self.basemap.urcrnrlon,self.dlon),color="0.5",
            linewidth=0.25, dashes=[1,1,0.1,1], labels=[0,0,0,1],fontsize=12)
        self.basemap.drawparallels(np.arange(self.basemap.llcrnrlat,
            self.basemap.urcrnrlat,self.dlat),color="0.5",
            linewidth=0.25, dashes=[1,1,0.1,1], labels=[1,0,0,0],fontsize=12)
        
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
        x,y = self.basemap(xrange, yrange)
        self.ax.set_xlim(x)
        self.ax.set_ylim(y)
        self.fig.canvas.draw()

    def pcolor(self, lon, lat, data, **kwargs):
        """
        pcolor field data onto our geographic plot
        
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
        dlon=lon*0;
        dlat=lat*0;
        dlon[:,0:-1]=lon[:,1:]-lon[:,0:-1]
        dlat[0:-1,:]=lat[1:,:]-lat[0:-1,:]
        x,y = self.basemap(lon-dlon*0.5,lat-dlat*0.5)
        self.pc = self.ax.pcolor(x,y,data,**kwargs)

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
        if label != None:
            self.cb.set_label(label)
        
        