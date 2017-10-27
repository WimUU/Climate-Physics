# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 12:04:03 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt
import netCDF4 as ncdf
import drifter_analysis as da

#%%


data = np.genfromtxt('01_prepped.csv', delimiter=',')
era = ncdf.Dataset("era_data.nc", mode = 'r')

#%%
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
plt.show()

ax= plt.subplot(111,polar = True)
plt.scatter(a.buoyangle, a.speed)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.grid(True)
ax.set_rlabel_position(180)
ax.tick_params(labelsize = 13)
plt.title('Wind Rose Plot of the speed of Buoy')
plt.savefig('windrose1.png',dpi=300)
plt.show()