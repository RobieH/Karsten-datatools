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
import os
from multiprocessing import Pool

def clearall():
    """clear all globals
    Not intended for user use.
    """
    for uniquevar in [var for var in globals().copy() if var[0] != "_" and var != 'clearall']:
        del globals()[uniquevar]
        
def energy_plot(data, show=True):
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
    if show:
        plt.show()
    
def el_plot(data, Map=False, show=True):
    """
    Plot the elevation for the region from the last time series
    
    @Parameters:
        data: the standard python data dictionary
        Map = {True, False} (optional): Optional argument.  If True,
            the elevation will be plotted on a map.  
    """
    trigrid = data['trigrid']
    plt.gca().set_aspect('equal')
    plt.tripcolor(trigrid, data['zeta'][-1,:])
    plt.colorbar()
    plt.title("Elevation")
    if Map:
        #we set the corners of where the map should show up
        llcrnrlon, urcrnrlon = plt.xlim()
        llcrnrlat, urcrnrlat = plt.ylim()
        #we construct the map.  Note that resolution serves to increase
        #or decrease the detail in the coastline.  Currently set to 
        #'i' for 'intermediate'
        m = Basemap(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, \
            resolution='i', suppress_ticks=False)
        #set color for continents.  Default is grey.
        m.fillcontinents(color='ForestGreen')
        m.drawmapboundary()
        m.drawcoastlines()
    if show:
        plt.show()

def speed_plot(data, Map=False, show=True):
    """
    Plot the speed in each element for the last time series
    
    @Parameters: 
        data: standard python data dictionary
        Map = {True, False} (optional): If True, the plot will include a map
        show = {True, False} (optional): If True, the plot will be 
            displayed.  If False, no plot will be displayed.
    """
    trigrid = data['trigrid']
    #see if speed  has already been calculated
    try:
        speed = data['speed']
    except KeyError:
        #calcuate speed
        sdata = calc_speed(data)
        speed = sdata['speed']
    #plot the speed for the last time series
    plt.gca().set_aspect('equal')
    plt.tripcolor(trigrid, speed[-1,:])
    plt.colorbar()
    if Map:
        #we set the corners of where the map should show up
        llcrnrlon, urcrnrlon = plt.xlim()
        llcrnrlat, urcrnrlat = plt.ylim()
        #we construct the map.  Note that resolution serves to increase
        #or decrease the detail in the coastline.  Currently set to 
        #'i' for 'intermediate'
        m = Basemap(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, \
            resolution='i', suppress_ticks=False)
        #set color for continents.  Default is grey.
        m.fillcontinents(color='ForestGreen')
        m.drawmapboundary()
        m.drawcoastlines()
    if show:
        plt.show()

def mass_tri_plot(data, savedir, name='plot', Type='speed', Map=False):
    """
    Plots all time series.
    """
    trigrid = data['trigrid']
    #get the data to plot
    try:
        toPlot = data[Type]
    except KeyError:
        print Type + " is not an element of data.  Please calculate it."
        raise Exception("Invalid dictionary entry")
    #set the variable as a global variable
    global plotvar 
    plotvar = toPlot
    global saveDir
    saveDir = savedir
    global grid
    grid = trigrid
    #see if the save directory exists, or make it
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    l = toPlot.shape[0]
    
    p = Pool(4)
    plt.gca().set_aspect('equal')
    p.map(save_plot, range(50))
    clearall()
    
def save_plot(i):       
        plt.tripcolor(grid, plotvar[i,:])
        plt.colorbar()
        plt.title('speed' + ' Time Series' + str(i))
        plt.savefig(saveDir + 'plot' + '_' + str(i) + '.png')
        plt.clf()