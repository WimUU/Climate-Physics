# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 21:47:25 2017

@author: Wim
"""

import numpy as np
import matplotlib.pyplot as plt

tag = '1'       # drifter number

coordinates    = np.load('coordinates'+tag+'.npy')
velocity       = np.load('velocity'+tag+'.npy')
angles_era_5   = np.load('angles_era_5'+tag+'.npy')
angles_era_int = np.load('angles_era_int'+tag+'.npy')
u_era_5        = np.load('u_era_5'+tag+'.npy')
v_era_5        = np.load('v_era_5'+tag+'.npy')
u_era_int      = np.load('u_era_int'+tag+'.npy')
v_era_int      = np.load('v_era_int'+tag+'.npy')
geostrophic    = np.load('geo_velocity_drifter'+tag+'.npy')

a_geostrophic = velocity-geostrophic[:231]

#%%

era_int_velocities = np.zeros((len(u_era_int),2))
era_5_velocities   = np.zeros((len(u_era_int),2))

for k in range(len(u_era_int)):
    era_int_velocities[k] = [u_era_int[k],v_era_int[k]]
    era_5_velocities[k]   = [u_era_5[k],v_era_5[k]]

#%%

full_angles_era_5   = np.zeros(len(u_era_int))
full_angles_era_int = np.zeros(len(u_era_int))

ageo_angles_era_5   = np.zeros(len(u_era_int))
ageo_angles_era_int = np.zeros(len(u_era_int))

for k in range(len(u_era_int)):
    v0 = era_int_velocities[k]                    # wind  
    v1 = velocity[k]                              # drifter
    angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    full_angles_era_int[k] = np.degrees(angle)    # Positive angle --> buoy movement to the right of the wind
    
for k in range(len(u_era_int)):
    v0 = era_5_velocities[k]                      # wind  
    v1 = velocity[k]                              # drifter
    angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    full_angles_era_5[k] = np.degrees(angle)      # Positive angle --> buoy movement to the right of the wind

for k in range(len(u_era_int)):
    v0 = era_int_velocities[k]                    # wind  
    v1 = a_geostrophic[k]                         # drifter
    angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    ageo_angles_era_int[k] = np.degrees(angle)    # Positive angle --> buoy movement to the right of the wind
    
for k in range(len(u_era_int)):
    v0 = era_5_velocities[k]                      # wind  
    v1 = a_geostrophic[k]                         # drifter
    angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    ageo_angles_era_5[k] = np.degrees(angle)      # Positive angle --> buoy movement to the right of the wind
    
#%%
    
r = np.corrcoef(velocity[:,0],geostrophic[:231][:,0])

test_r = np.corrcoef(velocity[:,0],u_era_int)

#%%
    
plt.figure(figsize=(12,8))
plt.title(r'Drifter_01 $u$ and measured $u_g$ (m/s) correlation',fontsize=25,y=1.02)
plt.plot(np.arange(231),velocity[:,0],linewidth=2.2,label='Drifter u-component (24h)')
plt.plot(np.arange(231),geostrophic[:231][:,0],linewidth=2.2,label='Geostrophic u-velocity')
plt.xlim((0,230))
plt.ylim((-1,1.75))
plt.grid(True)
plt.legend(fontsize=20)
plt.text(80,-0.75,'r = '+str(round(r[0,1],3)),color='red',fontsize=27.5,bbox=dict(facecolor='none', edgecolor='black', pad=10.0))
plt.xlabel('Day (starting from 2013-12-12)',fontsize=20)
plt.ylabel('Velocity (m/s)',fontsize=20)
plt.tick_params(labelsize=20)
