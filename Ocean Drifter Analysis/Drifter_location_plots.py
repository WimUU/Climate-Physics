# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:50:03 2017

@author: Wim
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

#%% SECTION A: Retrieve buoy GPS-locations

x_full = []
y_full = []

for m in [01,02,03,04,05,06]:
    data = genfromtxt('AAEdrifter_'+str(m)+'.csv', delimiter=',')
    
    x = np.zeros(len(data))
    y = np.zeros(len(data))
    
    for i in range(len(data)):
        x[i] = data[i,0] 
        y[i] = data[i,1]
    
    x = x.tolist()
    y = y.tolist() 
    x_full.append(x)
    y_full.append(y)
    
#%%
    
# create new figure, axes instances.
fig=plt.figure(figsize=(12.5,10))
ax=fig.add_axes([0.1,0.1,0.8,0.8])
m = Basemap(llcrnrlon=140.,llcrnrlat=-69,urcrnrlon=230.,urcrnrlat=-40.,\
            rsphere=(6378137.00,6356752.3142),\
            resolution='l',projection='merc')
# plot the trajectories as a line
m.plot(x,y,linestyle ='-',linewidth=2,latlon=True)
m.drawcoastlines()
m.fillcontinents()
# draw parallels
m.drawparallels(np.arange(-70,-40,10),labels=[1,1,0,1])
# draw meridians
m.drawmeridians(np.arange(0,360,30),labels=[1,1,0,1])
ax.set_title('Tracjetories of the 20 drifters released during the 2013-2014 AEE expedition.')
plt.show()


