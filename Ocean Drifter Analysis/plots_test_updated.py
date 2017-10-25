import numpy as np
import matplotlib.pylab as plt
from numpy import genfromtxt
import drifter_analysis_updated as da

#%%

data = genfromtxt('AAEdrifter_02 (copy).csv', delimiter=',')

# time_interval must be 6, 12, 24 or 48
a = da.drifter_analysis(data_set=data,time_interval=48,era_int=True,angle_calc=True)

#%%

b = da.drifter_analysis(data_set=data,time_interval=48,era_5=True,angle_calc=True)

#%%

test1 = b.angles_era_5
test2 = a.angles_era_int

#%%

diff = test1-test2
#
x = np.arange(len(test2))
plt.figure()
plt.plot(x,diff)

