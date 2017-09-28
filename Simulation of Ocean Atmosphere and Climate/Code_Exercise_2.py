# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:48:35 2017

@author: Wim
"""

import numpy as np
import matplotlib.pylab as plt

#%% SECTION A: Definitions of function at t = 0 (C_0)
#              Grid properties (length, time horizon, dx, dt)                
#              Wind-speed (u_0) 

L = 2500000                 # length of the domain (meters)
dx = 2500                   # stepsize (meters)
T = 250000                  # total time (seconds)
dt = 25                     # time stepsize (seconds)
u_0 = 10                    # m/s
time_steps  = T//dt         # number of iterations
grid_points = L//dx         # number of grid-points

# 'hat' function at t = 0
C_0 = np.zeros(grid_points)
for i in range(int(grid_points*1125/2500),int(grid_points*1375/2500)):
    C_0[i] = 1

#%% SECTION B: Definitions of finite differnce schemes with periodic BCs    

def euler_forward(dx,L,dt,T,u_0,C_0):
    # INPUT: Grid properties (dx,L,dt,T), velocity (u_0) and function at t=0 (C_0)
    # OUTPUT: C_0 at time T, evaluated by Euler Forward fin. difference scheme

    courant = u_0*dt/dx               
    C1 = C_0[:]                     # copy of input
    C2 = C1*0                       # list for t+1 time-step
    
    for i in range(T//dt):          # time-steps
        for j in range(L//dx):      # update for each grid-point
            C2[j] = C1[j] - courant*(C1[j] - C1[j-1])
        C1 = C2[:]                  # update C1 to t+1 time
        
    return C1  
    
def matsuno_upwind(dx,L,dt,T,u_0,C_0):
    # INPUT: Grid properties (dx,L,dt,T), velocity (u_0) and function at t=0 (C_0)
    # Matsuno scheme is a predictor-forward time stepping method.
    # Calculate t+1 values, then use t+1 in dCdx of fin. diff. scheme (EF)
    # OUTPUT: C_0 at time T, evaluated by Matsuno/EF fin. difference scheme

    courant = u_0*dt/dx   
    C1 = C_0[:]      
    C2 = C1*0                       # list for t+1 time-step
    matsuno = C_0*0                 # list for predictor-forward step
    
    for i in range(T//dt):          # time-steps
        for j in range(L//dx):      # predictor update for each grid-point
            matsuno[j] = C1[j] - courant*(C1[j] - C1[j-1])
        for j in range(L//dx):      # update for each grid-point
            C2[j] = C1[j] - courant*(matsuno[j] - matsuno[j-1])
        C1 = C2[:]                  # update C_0 to t+1 time
            
    return C1
    
def lax_wendroff(dx,L,dt,T,u_0,C_0):
    # INPUT: Grid properties (dx,L,dt,T), velocity (u_0) and function at t=0 (C_0)
    # OUTPUT: C_0 at time T, evaluated with lax_wendroff fin. diff. scheme

    coeff_1 = u_0*dt/dx/2         
    coeff_2 = (u_0*dt/dx)**2/2
    
    C1 = C_0[:]
    C2 = C1*0                      # list for t+1 time-step
    
    for i in range(T//dt):          # time-steps
        for j in range(L//dx):      # update for each grid-point
            if j == L//dx-1:        # periodic boundary condition
                C2[j] = (C1[j] - coeff_1*(C1[0]-C1[j-1])
                        + coeff_2*(C1[0]-2*C1[j]+C1[j-1]))
            else:
                C2[j] = (C1[j] - coeff_1*(C1[j+1]-C1[j-1])
                        + coeff_2*(C1[j+1]-2*C1[j]+C1[j-1]))
        C1 = C2[:]                  # update C_0 to t+1 time
            
    return C1
    
#%% SECTION B: Evaluation and plots of the three finite difference schemes

run_schemes = True  # set 'run_schemes = True' to initialize all fin.diff runs

if run_schemes:
    a = euler_forward(dx,L,dt,T,u_0,C_0)
    b = matsuno_upwind(dx,L,dt,T,u_0,C_0)
    c = lax_wendroff(dx,L,dt,T,u_0,C_0)

#%% SECTION C: Spectral method

def spectral_method(x,L,dt,T,u_0,C_0):
    # Function evaluates advection equation using the spectral method.
    # Spectral Fourier-coefficients are integrated in time, and the signal in
    # in time-domain is then reconstructed using said coefficients.
    # INPUT: Grid properties (dx,L,dt,T), velocity (u_0) and function at t=0 (C_0),
    # spatial grid onto which functions are projected (x = np.linspace(0,L,dx+1)).
    # OUTPUT: Advection of C_0 at time T
    
    coefficients = []       # to-be-filled list of two-valued lists containing coefficients
    J = len(C_0)

    for k in range(J//2+1): # create list of spectral fourier coefficients at t=0
        coefficients.append([alpha_k(C_0,k),beta_k(C_0,k)])
        
    for t in range(T//dt):  # time-step coefficients using Euler-Forward Matsuno
                            # alpha_0 and beta_0 are constant (skip those)  
        for k in range(1,J//2+1):
            alpha = coefficients[k][0]
            beta = coefficients[k][1]
            # Euler-Forward Matsuno scheme
            matsuno_alpha = alpha + dt*2*np.pi*u_0*k/L*beta
            matsuno_beta = beta - dt*2*np.pi*u_0*k/L*alpha
            # Update coefficients to time-step t=t+1
            coefficients[k][0] = alpha + dt*2*np.pi*u_0*k/L*matsuno_beta
            coefficients[k][1] = beta - dt*2*np.pi*u_0*k/L*matsuno_alpha
            
    y = C_0*0                 # (empty) domain onto which functions are projected
    for k in range(0,J//2+1): # Sum over all functions
        y = y + np.real((coefficients[k][0]
                        +coefficients[k][1]*1j)*np.exp(2*np.pi*k*x/L*1j))
    return y
    
def alpha_k(c,k):
    # INPUT: function (c), wave-number (k)
    # OUTPUT: alpha_k spectral fourier coefficient of c for wave-number k
    J = len(c)
    alpha = 0
    if k == 0:
        for j in range(0,J):
            alpha = alpha + c[j]
        alpha = 1/J*alpha
    else:
        for j in range(0,J):
            alpha = alpha + c[j]*np.cos(-2*np.pi*k*j/J)
        alpha = 2/J*alpha
    return alpha
    
def beta_k(c,k):
    # INPUT: function (c), wave-number (k)
    # OUTPUT: beta_k spectral fourier coefficient of c for wave-number k
    J = len(c)
    beta = 0 
    if k == 0:
        None
    else: 
        for j in range(0,J):
            beta = beta + c[j]*np.sin(-2*np.pi*j*k/J)
        beta = 2/J*beta
    return beta
    
#%% SECTION D: Evaluation of advection equation using Spectral Method
    
x = np.linspace(0,L,L//dx)
spectral = spectral_method(x,L,dt,T,u_0,C_0)

#%%

make_plot = True        # 'make_plot = True' to plot spectral method evaluation
save_plot = True
    
if make_plot:
    x = x*1e-3
    plt.figure(figsize=(15,9))
    plt.plot(x,a,label='Euler Forward',linewidth=3)
    plt.plot(x,b,label='Matsuno',linewidth=3)
    plt.plot(x,c,label='Lax-Wendroff',linewidth=3)
    plt.plot(x,spectral,label='Spectral Method',linewidth=3,color='black')
    plt.plot(x,C_0,label='C_0 (t=0)',linewidth=3,color = '0.65',linestyle='--')
    plt.title('Advection of C(x,t) evaluated at T = 2.5e5 s',fontsize=30)
    plt.legend(prop={'size': 25.0},loc=2)
    plt.xticks(size=25), plt.yticks(size=25)
    plt.xlabel('Distance (km)',fontsize=30), plt.ylabel('Function value',fontsize=30)
    plt.grid(True)
    if save_plot:
        plt.savefig('evaluation.eps', format='eps', dpi=500)
    