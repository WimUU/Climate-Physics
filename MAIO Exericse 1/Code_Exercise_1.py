# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 09:36:35 2017

@author: Wim
"""

import numpy as np
import scipy as sp
import matplotlib.pylab as plt
import statsmodels.tsa.api as smt
from statsmodels.graphics.tsaplots import plot_acf, plot_pacf

#%% SECTION A: Set-up AR2 data set

def AR2(a1=0,a2=0,noise=False,n=0,var=1.0,plot=False,acor_plot=False,out=False):
    
    # AR2: a1, a2 are the AR2 coefficients
    # noise=True introduces noise
    # n is the number of data points
    # var is the variance (=1 standard)
    # Optional: Plot of data, autocorrelation. Output of AR2 data set.
    
    x = np.random.normal(size=n)
    w = np.zeros(n)
    
    if noise:
        w = np.random.normal(size=n,scale=var)
    
    for i in range(n):
        x[i] = a1*x[i-1] + a2*x[i-2] + w[i]
        
    if plot:
        data_plot(x,a1,a2,var)
    if acor_plot:
        plot_acf(x), plot_pacf(x)
    if out:
        return x
        
def data_plot(data,a1,a2,var):
    
    # Function to plot AR2 data set
    
    dim = len(data)
    x = np.arange(dim)
    
    fig = plt.figure(figsize=(7,4))
    ax1 = fig.add_subplot(111)
    ax1.plot(x,data)
    ax1.tick_params(axis='both',which='major',labelsize=18)
    ax1.set_xlabel('Data point',fontsize=20)
    ax1.set_ylabel('Value',fontsize=20)
    ax1.set_title('AR(2). Var = '+str(var)+', a1 = '
                  +str(a1)+', a2 = '+str(a2),fontsize=20)
    ax1.xaxis.grid(True), ax1.yaxis.grid(True)

#%% SECTION B: Part 1 and 2 of exercise 1

# Function call creates and plots an AR2 data set.
# a1 = 0.75, a2 = 0, noise with variance 1 is introduced. 

# IMPORTANT: Autocorrelation output shows autocorrelation function as well as
#            95% confidence band (SHADED).

AR2(0.85,-0.25,noise=True,n=1000,var=25.0,plot=False,acor_plot=False)

#%% SECTION C
# Retrieve and prepare the data sets.
# Plot their respective autocorrelation and partial autocorrelation functions.

# Data array holds each subsequent data set in its entries

indices = ['a','b','c','d','e']
data = []

for i in indices:
    data.append(np.genfromtxt('data1'+i+'.txt'))
     
# Stability test:
# Set stability_check = True to plot figures.

stability_check = False
locs = [0,1,2,3,4]

if stability_check:
    for i in locs:
        
        # Moving average using convolution, autocovariance
        smoothed = np.convolve(data[i],np.ones(100)/100)
        x= np.arange(len(smoothed))
        plt.figure()
        plt.plot(x,smoothed)
        plt.xlabel('Data point')
        plt.title('Moving average data set ' +str(indices[i]))
        
        autocov = smt.acovf(data[i])
        x = np.arange(len(autocov))
        plt.figure()
        plt.plot(x,autocov)
    
# CONCLUSION: there was a linear trend in the fifth data set.
# Normalize and detrend (all) data

c = 0 
for d in data:
    # sp.signal.detrend removes linear trend from time series
    d = sp.signal.detrend(d)
    data[c] = (d - np.min(d))/(np.max(d) - np.min(d))
    c += 1
    
# Correlation_plots = True plots the auto- and partial autocorrelation plots

correlation_plots = False

if correlation_plots:
    for j in data:
        plot_acf(j,lags=40),plot_pacf(j,lags=40)

#%% SECTION E
# Fit the data to an ARIMA model:
# Data set indices [a,b,c,d,e] = [0,1,2,3,4], indices = locs (see prv. section).
# Set residual_plot = True to plot auto- and partial autocorrelation of residual.

# Fit data to ARIMA(p,d,q) model and retrieve the residual and coefficients:

arima_parameters = [(0,0,1),(0,0,2),(1,0,0),(2,0,0),(0,0,0)]
residual_plot = False
coefficients = []

for i in locs:    
    mdl = smt.ARIMA(data[i],order=arima_parameters[i]).fit()
    residual = mdl.resid
    # First entry of params is the mean of data set: skip it. 
    coefficients.append(mdl.params[1:])
    
    # No spikes in these plots: residual is white noise (desired).
    if residual_plot:
        plot_acf(residual,lags=40),plot_pacf(residual,lags=40)

#%% SECTION F
### OPTIONAL -- Brute force order (fitting) of the data -- OPTIONAL ###
# Set model index ('mdl_index') manually

brute_force = False

if brute_force:

    best_aic = np.inf 
    best_order = None
    best_mdl = None
    
    rng = range(3)
    biner = [0,1]
    mdl_index = 0
    
    for k in biner:
        for i in rng:
                for j in rng:
                    try:
                        tmp_mdl = smt.ARIMA(sp.signal.detrend(data[mdl_index]), order=(i,k,j)).fit(trend='nc')
                        tmp_aic = tmp_mdl.aic
                        if tmp_aic < best_aic:
                            best_aic = tmp_aic
                            best_order = (i,k,j)
                            best_mdl = tmp_mdl
                    except: continue
    
    # Print results, which includes parameter fits
    best_mdl.summary()  