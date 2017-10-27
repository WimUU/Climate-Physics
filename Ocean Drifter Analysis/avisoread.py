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

y=np.random.randn(10, 720, 1440)

latbounds = [ -60 , -50 ]
lonbounds = [ 150 , 180 ]
lats=np.linspace(-90,90,y.shape[1])
lons=np.linspace(0,360,y.shape[2])

# latitude lower and upper index
latli = np.argmin( np.abs( lats - latbounds[0] ) )
latui = np.argmin( np.abs( lats - latbounds[1] ) ) 

# longitude lower and upper index
lonli = np.argmin( np.abs( lons - lonbounds[0] ) )
lonui = np.argmin( np.abs( lons - lonbounds[1] ) )  

direcs = np.array([os.listdir('E:\Bryan\Downloads\AVISO')])
fig = plt.figure()
my_map=Basemap(llcrnrlon=lonbounds[0],llcrnrlat=latbounds[0],urcrnrlon=lonbounds[1],urcrnrlat=latbounds[1],resolution='l',projection ='merc')
my_map.drawcoastlines()
my_map.fillcontinents(color='gray')
my_map.drawparallels(np.arange(-80., -40., 10.), labels=[1,0,0,0], fontsize=10)
my_map.drawmeridians(np.arange(0., 360., 10.), labels=[0,0,0,1], fontsize=10)

lon,lat = np.meshgrid(lons[lonli:lonui],lats[latli:latui])
lons, lats = np.meshgrid(lons[lonli:lonui:4],lats[latli:latui:4])
cs = None
q = None
cbar = None

def init():
    return my_map                                                                                                                                                                                                                                                                                          

def animate(i):
    global cs
    global q
    global cbar
    if cs != None:
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
    gc.collect()
    print(i)
    return my_map
cs = None
q= None
cbar = None

#mywriter = animation.FFMpegWriter(codec="libx264",bitrate=-1)
anim = animation.FuncAnimation(fig, animate, init_func=init, frames=234, interval=100,save_count=0, repeat=False) #interval = number of milliseconds between frames
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
#my_map.quiver(lons,lats,u,v,scale=20,latlon=True)
#plt.show()

#cs = my_map.pcolormesh(lons,lats,u,cmap='seismic', shading='flat',latlon=True)
#cbar = my_map.colorbar(cs, location='bottom', pad="10%")
#cbar.set_label('$[ms^{-1}]$')
#plt.clim(-1,1)
