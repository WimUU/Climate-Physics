# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:54:00 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt
import random
from scipy.fftpack import fft, fftfreq, fftshift
           
#%% SECTION A: functions for initializing the desired time series.

def dirac_d(spacing, end_time):
    # INPUT: 'spacing' is years between time steps
    # 'end_time' length (in years) of data set
    x = np.zeros(int(end_time/spacing))
    x[len(x)/2] = 1
    return x

def sine_wave(wave_length,end_time):
    #INPUT: wave_length, length (yrs) of data set
    x = []
    
    for i in range(end_time):
        x.append(np.sin(2*np.pi*1/wave_length*i))
    return x

def sawtooth(spacing, end_time,noise=False,variations=False):
    #INPUT; spacing (average length saw-tooth), end_time length data set (yrs)
    # 'noise = True' introduces gaussian noise. 
    # 'variations=True' introduces +- 0.02myr variations in sawtooth periods

    
    x = []                          # empty list to which data is appended
    a = 0                           # counter for keeping track of variations
    b = [0,0,0,0,0,0,0,0,0,0]       # if 'variations=False'
    
    if variations:
        b = [20,20,20,20,20,-20,-20,-20,-20,-20]
        random.shuffle(b)           # randomly shuffle +- 0.02myr variations
    
    for i in range(int(end_time/spacing)):
        y = 1/(end_time+b[a])      # slope of sawtooth 
        for j in range(100+b[a]):   # variation-adjusted period length
            x.append(y*spacing*j) 
            
        a += 1
    if noise:                       # add noise to each data point
        for i in range(len(x)):     # adjust scale to increase nose amplitude
            x[i] = x[i] + np.random.normal(scale=10.0)
    return x
    
#%% SECTION B: Fast Fourier Transforms of sections A's time series.

q = sawtooth(100,1000,variations=True)

n = len(q)
T = 1
x = np.linspace(0.0,n*T,n)
yf = fft(q)
xf = np.linspace(0.0,1.0/(2.0*T),n//2)
plt.figure()
plt.plot(xf,1.0/n*abs(yf[0:n//2]))
plt.grid()
plt.show()