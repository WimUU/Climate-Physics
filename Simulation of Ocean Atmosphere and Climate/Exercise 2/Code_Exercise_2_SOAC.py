# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:48:35 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt

dx = 25         # kilometers
L = 2500        # kilometers
u_0 = 10        # m/s

#C_0 = np.zeros(100)
#
#for i in range(int(100*1125/2500),int(100*1375/2500)):
#    C_0[i] = 1

def euler_forward(dx,L,dt,T,u_0):
    # CFL-condition: abs(u_0*dt/dx <= 1)

    time_st = int(T/dt)         # number of iterations
    spatial = int(L/dx)         # number of grid-points
    courant = u_0*dt/dx         # courant number
    
    C_0 = np.zeros(spatial)
    for i in range(int(spatial*1125/2500),int(spatial*1375/2500)):
        C_0[i] = 1
                 
    C = np.zeros(spatial)       # list for next time step
    
    for i in range(time_st):
        for j in range(spatial):
            C[j] = C_0[j] - courant*(C_0[j] - C_0[j-1])
        C_0 = C[:]            
            
    return C_0
    
    
def matsuno_upwind(dx,L,dt,T,u_0):
    time_st = int(T/dt)         # number of iterations
    spatial = int(L/dx)         # number of grid-points
    courant = u_0*dt/dx         # courant number
    
    C_0 = np.zeros(spatial)
    for i in range(int(spatial*1125/2500),int(spatial*1375/2500)):
        C_0[i] = 1
    
    matsuno = np.zeros(spatial)
    C = np.zeros(spatial)
    
    for i in range(time_st):
        for j in range(spatial):
            matsuno[j] = C_0[j] - courant*(C_0[j] - C_0[j-1])
            C[j] = C_0[j] - courant*(matsuno[j] - C_0[j-1])
        C_0 = C[:]            
            
    return C_0
    
def lax_wendroff(dx,L,dt,T,u_0):
    time_st = int(T/dt)         # number of iterations
    spatial = int(L/dx)         # number of grid-points
    coeff_1 = u_0*dt/dx/2         
    coeff_2 = (u_0*dt/dx)**2/2
    
    C_0 = np.zeros(spatial)
    C = np.zeros(spatial)
    for i in range(int(spatial*1125/2500),int(spatial*1375/2500)):
        C_0[i] = 1
    
    for i in range(time_st):
        for j in range(spatial):
            if j == spatial-1:
                C[j] = (C_0[j] - coeff_1*(C_0[0]-C_0[j-1])
                        + coeff_2*(C_0[0]-2*C_0[j]+C_0[j-1]))
            else:
                C[j] = (C_0[j] - coeff_1*(C_0[j+1]-C_0[j-1])
                        + coeff_2*(C_0[j+1]-2*C_0[j]+C_0[j-1]))
        C_0 = C[:]            
            
    return C_0
    
def spectral_method():
    return None
    
a = matsuno_upwind(25000,2500000,250,1*250000,10)
b = euler_forward(25000,2500000,250,1*250000,10)
c = lax_wendroff(25000,2500000,250,1.0*250000,10)
plt.plot(np.arange(2500/25),a)
plt.plot(np.arange(2500/25),b)
plt.plot(np.arange(2500/25),c)
    


