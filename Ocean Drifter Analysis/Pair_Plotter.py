# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 14:50:03 2017

@author: Wim
"""
from __future__ import print_function

from mpl_toolkits.basemap import Basemap
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt

#%% SECTION A: Retrieve buoy GPS-locations

x_full = []     # latitude
y_full = []     # longitude

indices = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20']
path = 'C:/Users/Wim/Desktop/MAIO/MAIO_Drifter_project/'

for m in indices:
    data = genfromtxt(path+'AAEdrifter_'+m+'.csv', delimiter=',')
    
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
"""
    Record animation of buoy pair trajectories at 24 frames per second.
"""

plt.ioff()    # such that plots don't cram up iPython memory for displaying
q = 0
colors = ['black','red','green','blue','deeppink','orange','cyan','blueviolet','yellow','gray']

for i in [0,2,4,6,8,10,12,14,16,18]:       #[0,2,4,6,8,10,12,14,16,18]   
    ss1 = len(x_full[i])//(24.*3.)
    ss2 = len(x_full[i+1])//(24.*3.)

    for k in range(0,int(24.*3.)):        # ss = stepsize
    
        fig=plt.figure(figsize=(12,12))
        ax=fig.add_axes([0.1,0.1,0.8,0.8])
        m = Basemap(llcrnrlon=140.,llcrnrlat=-69,urcrnrlon=222.,urcrnrlat=-40.,\
        rsphere=(6378137.00,6356752.3142),\
        resolution='l',projection='merc')
        
        if i > 0:
            for j in [0,2,4,6,8,10,12,14,16,18][0:(i//2)]:
                    m.plot(x_full[j],y_full[j],linestyle ='-',linewidth=1.,color=colors[j//2],latlon=True)
                    m.plot(x_full[j+1],y_full[j+1],linestyle ='-',linewidth=1.,color=colors[j//2],latlon=True)
                    m.plot(x_full[j][len(x_full[j])-2:len(x_full[j])-1], y_full[j][len(x_full[j])-2:len(x_full[j])-1],latlon=True,marker="o",markersize=15,c=colors[j//2])
                    m.plot(x_full[j+1][len(x_full[j+1])-2:len(x_full[j+1])-1], y_full[j+1][len(x_full[j+1])-2:len(x_full[j+1])-1],latlon=True,marker="o",markersize=15,c=colors[j//2])  
        
        m.plot(x_full[i][:int((k+1)*ss1)],y_full[i][:int((k+1)*ss1)],linestyle ='-',linewidth=1.,color=colors[i//2],latlon=True)
        m.plot(x_full[i+1][:int((k+1)*ss2)],y_full[i+1][:int((k+1)*ss2)],linestyle ='-',linewidth=1.,color=colors[i//2],latlon=True)
        
        m.plot(x_full[i][int((k+1)*ss1)-2:int((k+1)*ss1)-1], y_full[i][int((k+1)*ss1)-2:int((k+1)*ss1)-1],latlon=True,marker="o",markersize=15,c=colors[i//2])
        m.plot(x_full[i+1][int((k+1)*ss2)-2:int((k+1)*ss2)-1], y_full[i+1][int((k+1)*ss2)-2:int((k+1)*ss2)-1],latlon=True,marker="o",markersize=15,c=colors[i//2])  
        
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

#%%

import os

os.system("ffmpeg -f image2 -r 24 -i C:/Users/Wim/Desktop/MAIO/MAIO_Drifter_project/Pair_Animation/buoy%5d.png -vcodec mpeg4 -y C:/Users/Wim/Desktop/MAIO/MAIO_Drifter_project/Pair_Animation/anim.mp4")
