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

direcs = np.array([os.listdir('E:\Bryan\Downloads\AVISO')])
fig = plt.figure()
my_map=Basemap(llcrnrlon=150,llcrnrlat=-70,urcrnrlon=220,urcrnrlat=-40,resolution='l',projection ='merc')
my_map.drawcoastlines()
my_map.fillcontinents(color='gray')
my_map.drawparallels(np.arange(-80., 90., 30.), labels=[1,0,0,0], fontsize=10)
my_map.drawmeridians(np.arange(0., 360., 30.), labels=[0,0,0,1], fontsize=10)

y=np.random.randn(10, 720, 1440) 

lats=np.linspace(-90,90,y.shape[1])
lons=np.linspace(0,360,y.shape[2])
lons, lats = np.meshgrid(lons,lats)
cs = None

def init():
    return my_map                                                                                                                                                                                                                                                                                          

def animate(i):
    global cs
    if cs != None:
        cs.remove()
    aviso = ncdf.Dataset('E:/Bryan/Downloads/AVISO/'+direcs[:,i][0], mode ='r')
    u = aviso.variables['u'][0,:,:]
    v = aviso.variables['v'][0,:,:]
    N = np.sqrt(u**2+v**2)
    cs = my_map.pcolormesh(lons,lats,u,cmap='seismic', shading='flat',latlon=True)
    cbar = my_map.colorbar(cs, location='bottom', pad="10%")
    cbar.set_label('$[ms^{-1}]$')
    plt.clim(-1,1)
#    cs = my_map.quiver(lons,lats,u/N,v/N,scale=15)
    gc.collect()
    print(i)
    return my_map
cs = None

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=10, interval=100,save_count=0) #interval = number of milliseconds between frames
anim.save('movie.mp4')

#aviso = ncdf.Dataset("dt_global_twosat_madt_uv_20131211_20140704.nc", mode = 'r')
##lats = aviso.variables['lat'][:]
##lons = aviso.variables['lon'][:]
#u = aviso.variables['u'][0,:,:]
#v = aviso.variables['v'][0,:,:]
#N=np.sqrt(u**2+v**2)
##lon,lat = np.meshgrid(lons,lats)
###xi,yi = my_map(lon,lat)
###
##
#m = Basemap(llcrnrlon=150,llcrnrlat=-70,urcrnrlon=220,urcrnrlat=-40,resolution='l')
#m.drawparallels(np.arange(-70., -45., 10.), labels=[1,0,0,0], fontsize=10)
#m.drawmeridians(np.arange(160., 260., 30.), labels=[0,0,0,1], fontsize=10)
#m.quiver(lons,lats,u,v,scale=15)
#plt.show()
###cs = m.pcolor(lon,lat,u,cmap='seismic')
#cs = m.pcolormesh(lon,lat,u,cmap='seismic')
#cbar = m.colorbar(cs, location='bottom', pad="10%")
#cbar.set_label('$[ms^{-1}]$')
#plt.clim(-1,1)
#m.drawcoastlines()
##m.fillcontinents(color='gray')
#plt.savefig('test.png',dpi=300)
#plt.show()
