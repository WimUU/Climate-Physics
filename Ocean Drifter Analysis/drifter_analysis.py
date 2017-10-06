# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:10:09 2017

@author: Wim
"""

import numpy as np
from geopy.distance import vincenty

class drifter_analysis(object):
    """
    Converts data from GPA-measurements to a list of velocities averaged over
    'time_interval' hours of buoy movement. 
    
    The imported data file should be structured as follows, per measurement/line: 
    # index 0: longitude
    #       1: latitude
    #       2: year
    #       3: month
    #       4: day
    #       5: hour
    """
    
    def __init__(self, data_set=None,time_interval=None,era_data=None,angle_calc=False):
        """
        If an input parameter is physical, use "physical" units, i.e. a diameter could be specified in meters.
        :param data_set: treated list of buoy GPS coordinates. List should start at
                         0 hours of the first day in consideration. Should only contain
                         complete daily measurements (i.e. list is a multiple of 24).
        :param time_interval: buoy data is recorded in hourly intervals. Velocity can be 
                         calculated by using different time intervals. For example,
                         time_interval = 12 calculates velocity from 12 hour interval changes
                         in location.
        """
        
        # intialize parameters
        self.data = data_set
        self.interval = time_interval
        
        # calculate coordinates and velocity data using function analysis()
        self.coordinates = None
        self.velocity = None
        self.speed = None
        self.analysis()    
        
        # ERA-interim data to be used in angle calculation
        self.era = era_data
        self.era_latitude = None
        self.era_longitude = None        
        self.era_wind = None
        self.angles = None
        if angle_calc:
            self.era_buoy_angle()
      
    def analysis(self):
        """
        Calculate velocities from input data. Returns list with consecutive
        velocities between period (hrs) of set interval.
        """
        
        velocities = np.zeros((len(self.data)//self.interval-1,2))
        speed = np.zeros(len(self.data)//self.interval-1)
        coordinates = np.zeros((len(self.data)//self.interval,2))
        
        for k in range(1,len(self.data)//self.interval):
            # u-velocity from difference in latitude/dt            
            # v-velocity from difference in longitude/dt
            x0 = self.data[(k-1)*self.interval][0]
            x1 = self.data[k*self.interval][0]
            y0 = self.data[(k-1)*self.interval][1]
            y1 = self.data[k*self.interval][1]
        
            # vincenty is a method to determine distance between points on earth
            u = vincenty((x0,y0),(x1,y0)).meters/(self.interval*60*60)
            v = vincenty((x0,y0),(x0,y1)).meters/(self.interval*60*60)

            
            if x1 < x0:
                u = -u   # Because vincenty always returns a positive number
            
            if y1 < y0:
                v = -v   # Because vincenty always returns a positive number
        
            velocities[k-1] = [u,v]
            coordinates[k-1] = [x0,y0]
            speed[k-1] = np.sqrt(u**2+v**2)
            
            # return last GPS points because velocity is 1 entry shorter
            # due to requiring two points to calculate velocity.
            if k == len(self.data)//self.interval-1:
                coordinates[k] = [x1,y1] 

        self.velocity = velocities
        self.coordinates = coordinates
        self.speed = speed
        
    def era_buoy_angle(self):
        """
        lorem ipsum
        """
        
        lat = self.era.variables['latitude'][:]
        lon = self.era.variables['longitude'][:]
        
        era_start = [2013.0, 12.0, 1.0]     # Date where ERA-data begins
        data_start = self.data[0][2:5]           # Date where 01_prepped begins
        
        skim_time_begin = int((data_start - era_start)[2]*4)  # 4 measurements per day ERA
        skim_time_end = int(skim_time_begin + len(self.data)/24*4) # 4 measurements per day ERA
        
        # Skim ERA-data to match time-interval of the measurements
        u_era = self.era.variables['u10'][skim_time_begin:skim_time_end,:,:]
        v_era = self.era.variables['v10'][skim_time_begin:skim_time_end,:,:]
                
        velocities_wind = np.zeros((len(self.data)//self.interval,2))

        for k in range(len(self.data)//self.interval):
            
            # buoy location relative to era-grid
            p = self.interval//6             # 6 because of ERA-interim spacing
            
            lon_test = self.coordinates[k][0]
            lat_test = self.coordinates[k][1]
            lon_diff = lon_test - 150
            lat_diff = -45 - lat_test
            
            lon_index = int(lon_diff//0.125) # Match buoy location with ERA-grid!
            lat_index = int(lat_diff//0.125) # Match buoy location with ERA-grid!
            
            # Average over four ERA-points because there are four per day
            
            u_culum = 0
            v_culum = 0   
            
            for m in range(p):            
                u_culum += u_era[int(k*p+m)][lat_index][lon_index]
                v_culum += v_era[int(k*p+m)][lat_index][lon_index]
                  
            u_culum = u_culum/p
            v_culum = v_culum/p

            velocities_wind[k] = (u_culum,v_culum)
            
        angles = np.zeros(len(self.velocity))
        
        for k in range(len(self.velocity)):
            v0 = velocities_wind[k]
            v1 = self.velocity[k]
            angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))    
            # Positive angle --> buoy movement to the right of the wind
            angles[k] = np.degrees(angle)
            
        self.era_wind = velocities_wind[:]
        self.era_latitude = lat[:]
        self.era_longitude = lon[:]
        self.angles = angles[:]