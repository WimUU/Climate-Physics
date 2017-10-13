import numpy as np
import matplotlib.pylab as plt
from numpy import genfromtxt
#import netCDF4 as ncdf
import drifter_test as da

#%%

data = genfromtxt('AAEdrifter_01 (copy).csv', delimiter=',')
#era = ncdf.Dataset("era_data.nc", mode = 'r')

# time_interval can be 12, 24 or 48 hours
a = da.drifter_analysis(data_set=data,time_interval=12)

#%%

#x = np.linspace(0,(len(a.angles)*time_interval//24),len(a.velocity))
#plt.figure(figsize=(12,8))
#plt.plot(x,a.angles,linewidth=2)
#plt.xlim(0,(len(a.angles)*time_interval//24))
#plt.title(r'Angle between wind and buoy velocity. Positive $\theta$ is to the right of the wind.',fontsize=17.5)
#plt.xlabel(r'Day (starting from 2013-12-12)',fontsize=17.5)
#plt.ylabel(r'$\theta$',fontsize=25)
#plt.grid()
#plt.show()
#
#ax = plt.subplot(111, polar = True)
#c = plt.scatter(a.angles, a.speed)
#ax.set_theta_zero_location('N')
#ax.set_theta_direction(-1)
#ax.grid(True)
#ax.set_rlabel_position(300)
#ax.tick_params(labelsize = 13)
#plt.show()
