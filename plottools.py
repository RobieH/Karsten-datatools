# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 21:20:34 2013

@author: robie

A number of functions for plotting various data features.
"""
import numpy as np
import matplotlib.tri as mplt
import matplotlib.pyplot as plt
from datatools import *
from mpl_toolkits.basemap import Basemap

def energy_plot(data):
    """
    Produces a plot of the kinetic energy as a function of time.  
    """
    #see if EK has already been calculated
    try:
        ek = data['ek']
    except KeyError:
        #calculate the energy as a function of time
        edata = calc_energy(data)
        ek = edata['ek']
    #plot the EK as a function of time.
    plt.plot(data['time'], ek)
    plt.title("Kinetic Energy vs. Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Kinetic Energy (J)")
    plt.show()
    
def el_plot(data, Map=False):
    trigrid = data['trigrid']
    plt.gca().set_aspect('equal')
    plt.tripcolor(trigrid, data['zeta'][-1,:])
    plt.colorbar()
    plt.title("Elevation")
    plt.xticks([data['lon'][0], data['lon'][-1]], ('cat', 'dog'))
    plt.yticks()   
    if Map:
        llcrnrlon, urcrnrlon = plt.xlim()
        llcrnrlat, urcrnrlat = plt.ylim()
        m = Basemap(llcrnrlon-3, llcrnrlat-3, urcrnrlon+3, urcrnrlat+3, \
            resolution='i')
        m.fillcontinents()
        m.drawmapboundary()
        m.drawcoastlines()
 
    plt.show()