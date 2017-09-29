# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 12:04:03 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt
from numpy import genfromtxt
import netCDF4 as ncdf
from geopy.distance import vincenty
import drifter_analysis as da

#%%

# index 0: longitude
#       1: latitude
#       2: year
#       3: month
#       4: day
#       5: hour

data = genfromtxt('01_prepped.csv', delimiter=',')

#%%
# Velocity measurements

interval = 24      # MUST BE MULTIPLE OF 3 (required for ERA-interim averaging)
velocities = []
coordinates = []

for k in range(1,len(data)//interval):
    
    x0 = data[(k-1)*interval][0]
    x1 = data[k*interval][0]
    y0 = data[(k-1)*interval][1]
    y1 = data[k*interval][1]

    u = vincenty((x0,y0),(x1,y0)).meters/(interval*60*60)
    v = vincenty((x0,y0),(x0,y1)).meters/(interval*60*60)
    
    if x1 < x0:
        u = -u
    
    if y1 < y0:
        v = -v

    velocities.append((u,v))
    coordinates.append((x0,y0))
    
    if k == len(data)//interval-1:
        coordinates.append((x1,y1))
        
#%%
        
print('The first measurement was taken on '+str(data[0][2])+'-'+str(data[0][3])+'-'+str(data[0][4]))
print('The data set ran for '+str(len(data)//24)+' consecutive days')
print('The velocity is calculated over '+str(interval)+' hour intervals.')
print('This results in '+str(len(velocities))+' velocity data points.')
    
#%%

a = da.drifter_analysis(data_set=data,time_interval=interval)

#%% Retrieve ERA-interim data

direc = "/Users/Wim/Desktop/MAIO/MAIO Drifter project/"
era = ncdf.Dataset(direc+"era_data.nc", mode = 'r')

lat = era.variables['latitude'][:]
lon = era.variables['longitude'][:]

#%%

era_start = [2013.0, 12.0, 1.0]     # Date where ERA-data begins
data_start = data[0][2:5]           # Date where 01_prepped begins

skim_time_begin = int((data_start - era_start)[2]*4)  # 4 measurements per day ERA
skim_time_end = int(skim_time_begin + len(data)/24*4) # 4 measurements per day ERA

# Skim ERA-data to match time-interval of the measurements
u_era = era.variables['u10'][skim_time_begin:skim_time_end,:,:]
v_era = era.variables['v10'][skim_time_begin:skim_time_end,:,:]

#%% Set up wind-velocities from ERA-data

velocities_wind = []

for k in range(len(data)//24):
    
    lon_test = data[k][0]
    lat_test = data[k][1]
    lon_diff = lon_test - 150
    lat_diff = -45 - lat_test
    
    lon_index = lon_diff//0.125 # Match buoy location with ERA-grid!
    lat_index = lat_diff//0.125 # Match buoy location with ERA-grid!
    
    # Average over four ERA-points because there are four per day
    u_meteo = (u_era[k][lat_index][lon_index] + u_era[1+k][lat_index][lon_index]
        + u_era[2+k][lat_index][lon_index] + u_era[4+k][lat_index][lon_index])/4
    v_meteo = (v_era[k][lat_index][lon_index] + v_era[1+k][lat_index][lon_index]
        + v_era[2+k][lat_index][lon_index] + v_era[4+k][lat_index][lon_index])/4
        
    velocities_wind.append((u_meteo,v_meteo))   # u_meteo and v_meteo = ERA
    
#%% We are interested in the angle between the buoy velocity and wind-velocity
    
angles = []

for k in range(len(velocities)):
    v0 = velocities_wind[k]
    v1 = velocities[k]
    angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))    
    # Positive angle --> buoy movement to the right of the wind
    angles.append(np.degrees(angle))
    
#%% Make plot
    
x = np.linspace(0,len(velocities)-1,len(velocities))
plt.figure(figsize=(12,8))
plt.plot(x,angles,linewidth=2)
plt.xlim(0,len(velocities))
plt.title(r'Angle between wind and buoy velocity. Positive $\theta$ is to the right of the wind.',fontsize=17.5)
plt.xlabel(r'Day (starting from 2013-12-12)',fontsize=17.5)
plt.ylabel(r'$\theta$',fontsize=25)
plt.grid()
    
#%% To, for example, retrieve GPS locations (see drifter_analysis class):
    
locations = a.coordinates
    
    
    
    
    
    
    
    