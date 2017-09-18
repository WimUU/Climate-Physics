# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:48:35 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt

dx = 2500         # kilometers
L = 2500000        # kilometers
dt = 25
T = 250000
u_0 = 10        # m/s
time_steps  = int(T/dt)         # number of iterations
grid_points = int(L/dx)         # number of grid-points

C_0 = np.zeros(grid_points)
for i in range(int(grid_points*1125/2500),int(grid_points*1375/2500)):
    C_0[i] = 1

def euler_forward(dx,L,dt,T,u_0):
    courant = u_0*dt/dx             # courant number                 
    time_steps  = int(T/dt)         # number of iterations
    grid_points = int(L/dx)         # number of grid-points
    C = np.zeros(grid_points)       # list for next time step

    C_0 = np.zeros(grid_points)
    for i in range(int(grid_points*1125/2500),int(grid_points*1375/2500)):
        C_0[i] = 1
    
    for i in range(time_steps):
        for j in range(grid_points):
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
        for j in range(spatial):
            C[j] = C_0[j] - courant*(matsuno[j] - matsuno[j-1])
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
    
#%% SECTION B: Plots of the three finite difference schemes
    
a = matsuno_upwind(2500,2500000,25,1.0*250000,10)
b = euler_forward(2500,2500000,25,1.0*250000,10)
c = lax_wendroff(2500,2500000,25,1.0*250000,10)
plt.plot(np.arange(25000/25),a)
plt.plot(np.arange(25000/25),b)
plt.plot(np.arange(25000/25),c)

#%% SECTION C: Spectral method

def spectral_method():
    
    return None

J = len(C_0)

def alpha_k(c,k):
    J = len(c)
    alpha = 0
    if k == 0:
        for j in range(0,J-1):
            alpha = alpha + c[j]
        alpha = 1/J*alpha
    else:
        for j in range(0,J-1):
            alpha = alpha + c[j]*np.cos(-2*np.pi*k*j/J)
        alpha = 2/J*alpha
    return alpha
    
def beta_k(c,k):
    J = len(c)
    beta = 0 
    if k == 0:
        None
    else: 
        for j in range(0,J-1):
            beta = beta + c[j]*np.sin(-2*np.pi*j*k/J)
        beta = 2/J*beta
    return beta
    
coefficients = []

for k in range(J//2):
    coefficients.append([alpha_k(C_0,k),beta_k(C_0,k)])
    
#%%
    
for t in range(T//dt):
    for k in range(1,J//2):
        alpha = coefficients[k][0]
        beta = coefficients[k][1]
        # time integration scheme
        matsuno_alpha = alpha + dt*2*np.pi*u_0*k/L*beta
        matsuno_beta = beta - dt*2*np.pi*u_0*k/L*alpha
        coefficients[k][0] = alpha + dt*2*np.pi*u_0*k/L*matsuno_beta
        coefficients[k][1] = beta - dt*2*np.pi*u_0*k/L*matsuno_alpha
    
#%%
    
x = np.linspace(0,L,dx+1)
y = x*0

for k in range(0,J//2):
    y = y + np.real((coefficients[k][0]+coefficients[k][1]*1j)*np.exp(2*np.pi*k*x/L*1j))
    
#%%
    
plt.plot(x,y)
    
