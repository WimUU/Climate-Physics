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
import os
import gc

folder = 'E:/Bryan/Downloads/'
direcs = np.array([os.listdir(folder + 'AVISO')])

coordinates_1 = np.load(folder + 'driftcoord/drifter_01.npy')
coordinates_2 = np.load(folder + 'driftcoord/drifter_02.npy')
velocity1b = np.load(folder + 'bijlagen/velocity1.npy')
velocityaviso = np.load(folder +'bijlagen/geo_velocity_drifter1.npy')
erav = np.load(folder + 'bijlagen/v_era_51.npy')
erau = np.load(folder + 'bijlagen/u_era_51.npy')

y=np.random.randn(10, 720, 1440)

latbounds = [ -61 , -49 ]
lonbounds = [ 150 , 181 ]
lats=np.linspace(-90,90,y.shape[1])
lons=np.linspace(0,360,y.shape[2])

# latitude lower and upper index
latli = np.argmin( np.abs( lats - latbounds[0] ) )
latui = np.argmin( np.abs( lats - latbounds[1] ) ) 

# longitude lower and upper index
lonli = np.argmin( np.abs( lons - lonbounds[0] ) )
lonui = np.argmin( np.abs( lons - lonbounds[1] ) )  

fig = plt.figure()
my_map=Basemap(llcrnrlon=lonbounds[0],llcrnrlat=latbounds[0]+0.3,urcrnrlon=lonbounds[1]-0.3,urcrnrlat=latbounds[1]-0.3,resolution='l',projection ='merc')
my_map.drawcoastlines()
my_map.fillcontinents(color='gray')
my_map.drawparallels(np.arange(-80., -40., 5.), labels=[1,0,0,0], fontsize=10)
my_map.drawmeridians(np.arange(0., 360., 10.), labels=[0,0,0,1], fontsize=10)


lon,lat = np.meshgrid(lons[lonli:lonui],lats[latli:latui])
lons, lats = np.meshgrid(lons[lonli:lonui:4],lats[latli:latui:4])
cs = None
q = None
cbar = None
x,y = my_map(0, 0)
punt1 = my_map.plot(x, y, 'ro', markersize=8)[0]
punt2 = my_map.plot(x, y, 'ro', markersize=8)[0]


def init():
    return my_map                                                                                                                                                                                                                                                                                          

def animate(i):
    global cs
    global q
    global cbar
#    global punt
    if cs != None:
        cbar.remove()
        cs.remove()
        q.remove()
    aviso = ncdf.Dataset('E:/Bryan/Downloads/AVISO/'+direcs[:,i][0], mode ='r')
    u = aviso.variables['u'][0,latli:latui,lonli:lonui]
    v = aviso.variables['v'][0,latli:latui,lonli:lonui]
    us = u[::4,::4]
    vs = v[::4,::4]
    N = np.sqrt(u**2+v**2)
    cs = my_map.pcolormesh(lon,lat,N,cmap='jet', shading='flat',latlon=True)
    cbar = my_map.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label('$[ms^{-1}]$')
    plt.clim(0,1.5)
    q = my_map.quiver(lons,lats,us,vs, scale =15, latlon=True)
    x,y = my_map(coordinates_1[i,0],coordinates_1[i,1])
    punt1.set_data(x,y)
    x2,y2 = my_map(coordinates_2[i,0],coordinates_2[i,1])
    punt2.set_data(x2,y2)
    gc.collect()
    print(i)
    return my_map
cs = None
q= None
cbar = None


# kies aantal frames niet langer dan 232!
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=120, interval=100,save_count=0, repeat=False) #interval = number of milliseconds between frames
#anim.save('movie.mp4', bitrate=500,dpi=300)
anim.save('filmpje2.mp4',bitrate=500, dpi=300)

#%%
#aviso = ncdf.Dataset('E:/Bryan/Downloads/AVISO/'+direcs[:,1][0], mode ='r')
#fig = plt.figure()
#my_map=Basemap(llcrnrlon=150,llcrnrlat=-70,urcrnrlon=220,urcrnrlat=-40,resolution='l',projection ='merc')
#my_map.drawcoastlines()
#my_map.fillcontinents(color='gray')
#my_map.drawparallels(np.arange(-80., 90., 30.), labels=[1,0,0,0], fontsize=10)
#my_map.drawmeridians(np.arange(0., 360., 30.), labels=[0,0,0,1], fontsize=10)
#u = aviso.variables['u'][0,:,:]
#v = aviso.variables['v'][0,:,:]
#N = np.sqrt(u**2+v**2)
#x,y = my_map(160,-55)
#my_map.plot(x,y,marker='D', color = 'm')
#plt.show()

#cs = my_map.pcolormesh(lons,lats,u,cmap='seismic', shading='flat',latlon=True)
#cbar = my_map.colorbar(cs, location='bottom', pad="10%")
#cbar.set_label('$[ms^{-1}]$')
#plt.clim(-1,1)

#%%

np.corrcoef(velocity1b[:,1],velocityaviso[:231,1])

plt.plot(velocity1b[:,1])
plt.plot(velocityaviso[:,1])
plt.grid()
#plt.plot(erav)
plt.show()

np.corrcoef(velocity1b[:,1]-velocityaviso[:231,1],erav)
np.corrcoef(velocity1b[:,0]-velocityaviso[:231,0],erau)
