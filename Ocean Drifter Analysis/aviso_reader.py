# -*- coding: utf-8 -*-
"""
Created on Sun Oct  8 10:38:10 2017
@author: Bryan
"""

import numpy as np
import matplotlib.pylab as plt
from matplotlib import animation
import netCDF4 as ncdf
from mpl_toolkits.basemap import Basemap
import copy as copy
import os
import gc
     
#%%

folder = 'C:/Users/Wim/Desktop/maio/maiodrifterproject/'
direcs = np.array([os.listdir(folder+'aviso/')])

coordinates_1 = np.load(folder+'drifter_01.npy')
coordinates_2 = np.load(folder+'drifter_02.npy')

#%% Initialize grid indices

data = ncdf.Dataset(folder+'aviso/'+direcs[0,0], mode ='r')
lat = data.variables['lat']
lon = data.variables['lon']

lat_1 = -40
lat_2 = -70
lon_1 = 150
lon_2 = 220

lat_end = 0
for j in lat:
    if j > lat_1:
        break
    else:
        lat_end += 1

lat_start = 0
for j in lat:
    if j > lat_2:
        break
    else:
        lat_start += 1  
        
lon_end = 0
for j in lon:
    if j > lon_2:
        break
    else:
        lon_end += 1

lon_start = 0
for j in lon:
    if j > lon_1:
        break
    else:
        lon_start += 1 
#%%

lat = data.variables['lat'][lat_start:lat_end]        
lon = data.variables['lon'][lon_start:lon_end]

#%%

u_field = np.zeros((lat_end-lat_start,lon_end-lon_start,len(direcs[0,:])))
v_field = np.zeros((lat_end-lat_start,lon_end-lon_start,len(direcs[0,:])))

#%%

for j in range(len(direcs[0,:])):
    
    data = ncdf.Dataset(folder+'aviso/'+direcs[0,j], mode ='r')   
    u_field[:,:,j] = data.variables['u'][0,lat_start:lat_end,lon_start:lon_end]
    v_field[:,:,j] = data.variables['v'][0,lat_start:lat_end,lon_start:lon_end]

#%%

velocity_indices = np.zeros((len(coordinates_2),2))

for k in range(len(coordinates_2)):
    count_lat = 0
    count_lon = 0
    
    for j in lon:
        if j > coordinates_2[k][0]:
            break
        else:
            count_lon += 1    
    
    for j in lat:
        if j > coordinates_2[k][1]:
            break
        else:
            count_lat += 1
            
    if np.abs(lat[count_lat+1] - coordinates_2[k][1]) < np.abs(lat[count_lat] - coordinates_2[k][1]):
        count_lat += 1
    if np.abs(lon[count_lon+1] - coordinates_2[k][0]) < np.abs(lon[count_lon] - coordinates_2[k][0]):
        count_lon += 1
        
    velocity_indices[k] = [count_lat,count_lon]
    
#%%
    
save_geo_velocity = False
    
if save_geo_velocity:
    geo_velocity = np.zeros((len(coordinates_2),2))
    
    # coordinates_2 starts on the beginning of the 12th. Hence the k+1 in the u,v-field arrays
    for k in range(len(coordinates_2)):
        geo_velocity[k] = [u_field[velocity_indices[k][0],velocity_indices[k][1],k+1],v_field[velocity_indices[k][0],velocity_indices[k][1],k+1]]
    
    np.save(folder+'geo_velocity_drifter2.npy',geo_velocity)

#%%

q = 0
for k in range(len(u_field[0,0,:])):        # len(u_field[0,0,:])
    
    fig=plt.figure(figsize=(12,12))
    ax=fig.add_axes([0.1,0.1,0.8,0.8])
    m = Basemap(llcrnrlon=151.,llcrnrlat=-69,urcrnrlon=219.,urcrnrlat=-41.,\
    rsphere=(6378137.00,6356752.3142),\
    resolution='l',projection='merc')
                  
    temp_field = copy.copy(u_field[::-1,:,k])  
    temp_field = np.ma.masked_where((temp_field) < -100,temp_field)
    m.pcolormesh(lon,lat,temp_field,latlon=True,vmin=-0.5,vmax=0.5,cmap=plt.cm.bwr)  
        
    m.drawcoastlines()
    m.fillcontinents()
    m.drawparallels(np.arange(-70,-40,10),labels=[1,1,0,1])
    m.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1])
    ax.set_title('Trajectories of 10 drifter pairs released during the 2013-2014 AEE expedition.',fontsize=17.5)
        
    marker_c = 5 - len(str(q))
    marker = marker_c*'0'+str(q)
    plt.savefig(('buoy'+marker+'.png'),dpi=300)
        
    plt.close()
    q += 1

os.system("ffmpeg -f image2 -r 10 -i C:/Users/Wim/Desktop/maio/maiodrifterproject/buoy%5d.png -s 1920x1080 -vcodec libx264 -pix_fmt yuv420p -y C:/Users/Wim/Desktop/maio/maiodrifterproject/geostrophic.mp4")