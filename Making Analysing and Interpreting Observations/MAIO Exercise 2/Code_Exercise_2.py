# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:54:00 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt
import random

#%% SECTION A: Fourier transform algorithm
    
def fourier_transform(x):
    # INPUT: Array (1-Dimensional)
    # OUTPUT: Non-scaled array with spectral analysis of input (time) series
    x = np.asarray(x, dtype=float)  # convert to array (computationally beneficial)
    N = len(x)
    n = np.arange(N)                # harmonics
    k = n.reshape((N, 1))           # wave-number
    M = np.exp(-2j*np.pi*k*n/N)     # for each k evaluate for all n
    return np.dot(M, x) 
           
#%% SECTION B: functions for initializing the desired time series.

def unit_impulse(spacing, length,make_plot=False):
    # INPUT: 'spacing' is kyrs between time steps
    # length (in kyrs) of data set
    # OUTPUT: unit impulse halfway of the length of the time series
    x = np.zeros(int(length/spacing))
    x[len(x)//2] = 1
    
    if make_plot:
        n = len(x)
        yf = fourier_transform(x)
        xf = np.linspace(0.0,1.0/(2.0*spacing),n//2)
        plt.figure()
        plt.plot(xf,1.0/n*abs(yf[0:n//2]))
        plt.xlim((0,0.05))
        plt.title('Frequency domain',fontsize=15)
        plt.ylabel('DFT Value',fontsize=15), plt.xlabel('Frequency (per kyrs)',fontsize=15)
        plt.grid()
        plt.show()
        
    return x

def sine_wave(period,length,spacing, make_plot=False):
    #INPUT: period (kyrs), length (kyrs), (sample) spacing (kyrs)
    # 'int_pwr = True' plots the integration power function in the frequency domain
    #OUTPUT: sine-wave with amplitude 1 with its assigned period and length
    x = np.arange(length)
    x = np.sin(2*np.pi*1/period*x)

    if make_plot:
        n = len(x)
        yf = fourier_transform(x)
        xf = np.linspace(0.0,1/(2.0*spacing),n//2)
        plt.figure()
        plt.plot(xf,1.0/n*abs(yf[0:n//2]))
        plt.xlim((0,0.1))
        plt.title('Frequency domain',fontsize=15)
        plt.ylabel('DFT Value',fontsize=15), plt.xlabel('Frequency (per kyrs)',fontsize=15) 
        plt.grid()
        plt.show()
        
    return x

def sawtooth(period,length,spacing,noise=False,variations=False,make_plot=False):
    #INPUT: period (average, kyrs) end_time length data set (yrs)
    # 'noise = True' introduces gaussian noise. 
    # 'variations=True' introduces +- 0.02myr variations in sawtooth periods
    # 'int_pwr = True' plots the integration power function in the frequency domain
    # OUTPUT: Sawtooth function with amplitude 1

    x = []                          # empty list to which data is appended
    a = 0                           # counter for keeping track of variations
    b = [0,0,0,0,0,0,0,0,0,0]       # if 'variations=False'
    
    if variations:
        b = [20,20,20,20,20,-20,-20,-20,-20,-20]
        random.shuffle(b)           # randomly shuffle +- 0.02myr variations
    
    for i in range(int(length/period)):
        y = 1/(period+b[a])      # slope of sawtooth 
        for j in range(100+b[a]):   # variation-adjusted period length
            x.append(y*j) 
        a += 1
        
    if noise:                       # add noise to each data point
        for i in range(len(x)):     # adjust scale to increase nose amplitude
            x[i] = x[i] + np.random.normal(scale=.1)
            
    if make_plot:
        n = len(x)
        yf = fourier_transform(x)
        xf = np.linspace(0.0,1.0/(2.0*spacing),n//2)
        plt.figure()
        plt.plot(xf,1.0/n*abs(yf[0:n//2]))
        plt.title('Frequency domain',fontsize=15)
        plt.ylabel('DFT Value',fontsize=15), plt.xlabel('Frequency (per kyrs)',fontsize=15)
        plt.grid()
        plt.show()
        
    return x

def signal_plotter(s):
    x = np.arange(len(s))
    plt.figure()
    plt.plot(x,s)
    plt.title('Time domain',fontsize=15)
    plt.ylabel('f(t)',fontsize=15), plt.xlabel('t (kyrs)',fontsize=15)
    plt.grid()
    plt.show()
    
#%% SECTION C: Fourier Transforms of sections A's time series.
    
signal_plot = False # Set 'signal_plot = True' for time-domain plots

# Question 1. 0.6 Myr time-series (600 kyrs). impulse halfway at 300 kyrs.

impulse = unit_impulse(1,600,make_plot=False)

# Question 2: 1000 kyrs sine-wave with wavelength 100 kyrs, amplitude 1

sinwav = sine_wave(100,1000,1,make_plot=False)

# Question 3: 1000 kyrs sawtooth function with wavelength 100 kyrs, amplitude 1

saw = sawtooth(100,1000,1,make_plot=False)

int_pwr = True

if int_pwr:
    sin_transform = fourier_transform(sinwav)
    n = len(sin_transform)
    int_pwr_sin = np.sum(1.0/n*abs(sin_transform))

    saw_transform = fourier_transform(saw)
    n = len(saw_transform)
    int_pwr_saw = np.sum(1.0/n*abs(saw_transform))

# Question 5: 1000 kyrs sawtooth time series with randomly 20 kyrs variations 
# in period. Amplitude 1. 

saw_var = sawtooth(100,1000,1,variations=True,make_plot=False)

# Question 6: Repeat question 4 but with significant noise added to ampltiude

saw_var2 = sawtooth(100,1000,1,noise=True,variations=True,make_plot=False)

if signal_plot:
    signal_plotter(impulse)
    signal_plotter(sinwav)
    signal_plotter(saw)
    signal_plotter(saw_var)
    signal_plotter(saw_var2)
