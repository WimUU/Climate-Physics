# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 12:04:03 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt
from numpy import genfromtxt
import netCDF4 as ncdf
import drifter_analysis as da

#%%

#Even testen of GitHub lekker werkt.

data = genfromtxt('01_prepped.csv', delimiter=',')
era = ncdf.Dataset("/Users/Wim/Desktop/MAIO/MAIO Drifter project/era_data.nc", mode = 'r')

time_interval = 48
a = da.drifter_analysis(data_set=data,time_interval=time_interval,era_data = era, angle_calc=True)

#%%

x = np.linspace(0,(len(a.angles)*time_interval//24),len(a.velocity))
plt.figure(figsize=(12,8))
plt.plot(x,a.angles,linewidth=2)
plt.xlim(0,(len(a.angles)*time_interval//24))
plt.title(r'Angle between wind and buoy velocity. Positive $\theta$ is to the right of the wind.',fontsize=17.5)
plt.xlabel(r'Day (starting from 2013-12-12)',fontsize=17.5)
plt.ylabel(r'$\theta$',fontsize=25)
plt.grid()
