# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 17:52:26 2017

@author: Wim
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Taylor import TaylorDiagram 

#%% SECTION A: Global parameters

latitude = 52                              # latitudonal coordinate coast (deg)
omega = 7.2792e-5                          # angular velocity Earth (s^-1)
f = 2*omega*np.sin(latitude*np.pi/180)     # coriolis parameter (s^-1)
rho = 1.25                                 # density air (kg/m3)

#%% SECTION B: Retrieve data and construct list with geostrophic wind

measurements = pd.read_excel('data2.xlsx')
data = []

for index in measurements:
    data.append(index)
    
c = 0
for index in data:
    data[c] = measurements[index]
    c += 1
    
# Construct list with geostrophic winds using geostrophic balance. 
# Geostrophic balance:
# u = -1/(f*rho)*p_y, v = 1/(f*rho)*p_x
    
u_g = -1/(f*rho)*data[3]/1000              # correct Pa/km to Pa/m
v_g = 1/(f*rho)*data[2]/1000               # correct Pa/km to Pa/m

# 'mean_hours = int()' is the number of hrs over which geostr. wind is averaged
p = len(data[3])
mean_hours = int(24)

for j in range(int(p/mean_hours)):
    for i in range(mean_hours):
        u_g[i+int(mean_hours*j)]=( 
               -1/(f*rho)*np.mean(data[3][mean_hours*j:mean_hours*(j+1)])/1000)
        v_g[i+int(mean_hours*j)]=( 
                1/(f*rho)*np.mean(data[2][mean_hours*j:mean_hours*(j+1)])/1000)

u_g = u_g.tolist()                         # from pandas dataframe to list
v_g = v_g.tolist()                         # from pandas dataframe to list


# Alternative geostrophic wind: Daily mean (as used in the report)
u_g2 = []
v_g2 = []

for k in range(2):
    for i in range(3600//360*24):
        u_g2.append(u_g[24*k])
        v_g2.append(v_g[24*k])
        
        
u_g = u_g2[:]
v_g = v_g2[:]    

#%% SECTION C: Analytical solution and finite difference scheme functions

def analytical_solution(A,phi,dt,T):
    
    time = np.linspace(0,T*3600,T*3600/dt)    
    C1 = - A*f/(rho*(f**2 - omega**2))
    C2 = 0                                  # boundary condition u_0 = 0
    C3 = A*omega/(rho*(f**2 - omega**2))
    
    return C1*np.sin(f*time) + C2*np.cos(time) + C3*np.sin(omega*time + phi)
    
def num_model(A,phi,dt,T,u_g,v_g,rd_constant,init_cond=False,geostr=False):
    # INPUT: 
    # A: 
    # phi:
    # Total time T (hrs), stepsize dt (s)
    # Geostrophic wind u_g, v_g
    # Rayleigh damping coefficient, rd_constant
    # 'init_cond = True' for initial condition (velocities at t=0)
    # 'geostr = True' for inclusion of geostrophic wind
    # OUTPUT:
    # lists of seabreeze velocity field u and v
    # List of temporal domain (time)

    n_steps = (T*3600)//dt                  # number of time steps 
    time = np.linspace(0,T*3600,T*3600/dt)  # time discretization
    u = np.zeros(n_steps)                   # empty list for u (m/s)
    v = np.zeros(n_steps)                   # empty list for v (m/s)
      
    if init_cond:                           # initial velocities, 0 otherwise
        u[0] = -8
        v[0] = -7.5

    # set up interpolated list of geostrophic wind
    # list of geostrophic values has length 48. It must be expanded to fit
    # the dimensions of discretized time (length T*3600/dt)

    u_g = u_g[:]                            # interal copy of geowind-lists
    v_g = v_g[:]
    
    u_geo = u_g[:]
    v_geo = v_g[:]
     
    # Implementation of Runge-Kutta 4 (RK4) scheme. 
    # Set up arrays for RK4
    u_1 = u*0
    u_2 = u*0
    u_3 = u*0
    u_4 = u*0
    v_1 = v*0
    v_2 = v*0
    v_3 = v*0
    v_4 = v*0
                
    for t in range(n_steps-1):              # minus one step due to 0th index
        # Basic RK4 recepy for equation of the form u_t = f(u)
        u_1[t] = u[t] + dt/2*(f*(v[t]-v_geo[t]) 
                      - A/rho*np.cos(omega*time[t] + phi) - rd_constant*u[t])
        v_1[t] = v[t] - dt/2*(f*(u[t]-u_geo[t]) + rd_constant*v[t])
        
        u_2[t] = u[t] + dt/2*(f*(v_1[t]-v_geo[t]) 
                      - A/rho*np.cos(omega*time[t] + phi) - rd_constant*u_1[t])
        v_2[t] = v[t] - dt/2*(f*(u_1[t]-u_geo[t]) + rd_constant*v_1[t])
        
        u_3[t] = u[t] + dt*(f*(v_2[t]-v_geo[t]) 
                      - A/rho*np.cos(omega*time[t] + phi) - rd_constant*u_2[t])
        v_3[t] = v[t] - dt*(f*(u_2[t]-u_geo[t]) + rd_constant*v_2[t])
        
        u_4[t] = u[t] - dt/2*(f*(v_3[t]-v_geo[t]) 
                      - A/rho*np.cos(omega*time[t] + phi) - rd_constant*u_3[t])
        v_4[t] = v[t] + dt/2*(f*(u_3[t]-u_geo[t]) + rd_constant*v_3[t])
                      
        u[t+1] = (u_1[t]+2*u_2[t]+u_3[t]-u_4[t])/3 
        v[t+1] = (v_1[t]+2*v_2[t]+v_3[t]-v_4[t])/3
        
    return [u,v,time]                       # u,v evaluated at time T

#%% SECTION D: Calibration of 'A' and 'phi' parameters

def estimate(A,phi,t):
    # INPUT: A and phi (see above), and t (hrs)
    # OUTPUT: r.h.s. of equation 1.109 from section 1.17 in Pa/m
    return A*np.cos(omega*t*3600 + phi)

daily_average_plot = False                # set to 'True' to produce plot

if daily_average_plot:
    # line phase and amplitude of the two signals up to 'calibrate'
    # l.h.s. and r.h.s. of equation 1.109
    time = np.arange(24)                  # twenty four hours 
    plt.figure(figsize=(10,6))
    plt.plot(time, data[6][24:]/1000,label='Measured two-day average',linewidth=3)
    plt.plot(time,estimate(0.0013,5.2,time),label=r'Eq. 1.106 (AD) for A = 1e-3 Pa/m, $\phi$ = 5.5 Rad',linewidth=3,color='red')
    plt.legend(prop={'size': 17.5},loc=3)
    plt.xlim(0,23)
    plt.xlabel('Time (hours)',fontsize=22) 
    plt.ylabel(r'$\partial p / \partial x$ (Pa/m)',fontsize=22) 
    plt.grid()
    plt.xticks(fontsize=20), plt.yticks(fontsize=15)
    plt.title(r'Two-day average $\partial p / \partial x$ vs. calculated $\partial p / \partial x$',fontsize=22)
#    plt.savefig('fig2.eps', format='eps', dpi=500)

#%% SECTION F: Output 
    
a = num_model(0.0013,5.2,360,48,u_g,v_g,0.00016,init_cond=True,geostr=True)
#u_g = np.zeros(len(u_g))
#v_g = np.zeros(len(v_g))
#a = num_model(0.001,0,360,48,u_g,v_g,.0,init_cond=False,geostr=False)

make_plot = True

if make_plot: 
    inter_data_u = []
    n_points = len(a[2])
    interpolate = n_points//len(data[4])
    for i in range(len(data[4])):
        for j in range(interpolate):
            if j == 0:
                inter_data_u.append(data[4][i])
            else:
                inter_data_u.append(None)
                
    fig = plt.figure(figsize=(10,6))
    ax1 = fig.add_subplot(111)
    ax1.plot(a[2]/3600,a[0],label='Model (RK4)',linewidth=2)#,color='black')
    plt.scatter(a[2]/3600,inter_data_u,color='black',label='Measurements')
#    ax1.plot(a[2]/3600,analytical_solution(0.001,0,360,48),label='Analytical solution',linewidth=3,color='red',linestyle='--')
    ax1.tick_params(axis='both',which='major',labelsize=18)
    ax1.set_xlabel('Time (hrs)',fontsize=20)
    ax1.set_ylabel('u (m/s)',fontsize=20)
    ax1.set_title('Idealized seabreaze at IJmuiden',fontsize=20)
    ax1.set_title('Modelled and measured seabreeze at IJmuiden',fontsize=20)
    ax1.set_xlim((0,48))
    plt.legend(prop={'size': 17.5},loc=2)
    ax1.xaxis.grid(True), ax1.yaxis.grid(True)
#    plt.savefig('seabreeze.eps', format='eps', dpi=500)
    
#%% Taylor Diagram

# Reference std
stdrefs = dict(model=np.std(data[4]))

# Sample std,rho: Be sure to check order and that correct numbers are placed!
samples = dict(model=[ [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00008,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00008,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 0.8e-4"],
                       [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00012,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00012,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 1.2e-4"],
                       [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00016,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00016,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 1.6e-4"],
                       [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00020,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00020,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 2.0e-4"],
                       [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00024,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00024,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 2.4e-4"],
                       [np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00028,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00028,init_cond=True,geostr=True)[0][::10])[1,0], "$\lambda$ = 2.8e-4"]])

# Colormap (see http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps)
colors = plt.matplotlib.cm.Set1(np.linspace(0,1,len(samples['model'])))

# Here set placement of the points marking 95th and 99th significance
# levels. For more than 102 samples (degrees freedom > 100), critical
# correlation levels are 0.195 and 0.254 for 95th and 99th
# significance levels respectively. Set these by eyeball using the
# standard deviation x and y axis.

rects = dict(model=221)

fig = plt.figure(figsize=(6,5))

for season in ['model']:

    dia = TaylorDiagram(stdrefs[season], fig=fig,
                        label='Measurements')

    # Add samples to Taylor diagram
    for i,(stddev,corrcoef,name) in enumerate(samples[season]):
        dia.add_sample(stddev, corrcoef,
                       marker='$%d$' % (i+1), ms=13, ls='',
                       #mfc='k', mec='k', # B&W
                       mfc=colors[i], mec=colors[i], # Colors
                       label=name)

    # Add RMS contours, and label them
    contours = dia.add_contours(levels=5, colors='0.5') # 5 levels
    dia.ax.clabel(contours, inline=1, fontsize=12, fmt='%.1f')
    # Tricky: ax is the polar ax (used for plots), _ax is the
    # container (used for layout)
    dia._ax.set_title('Model performance',fontsize=15)

# Add a figure legend and title. For loc option, place x,y tuple inside [ ].
# Can also use special options here:
# http://matplotlib.sourceforge.net/users/legend_guide.html

fig.legend(dia.samplePoints,
           [p.get_label() for p in dia.samplePoints],
           numpoints=1, prop=dict(size=10))

fig.tight_layout()

#plt.savefig('test_taylor_4panel.eps')

#%% Model and data meta-data

print(np.std(num_model(0.0013,5.2,360,48,u_g,v_g,0.00016,init_cond=True,geostr=True)[0][::10]), np.corrcoef(data[4],num_model(0.0013,5.2,360,48,u_g,v_g,0.00016,init_cond=True,geostr=True)[0][::10])[1,0])