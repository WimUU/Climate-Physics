# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 15:10:09 2017

@author: Wim
"""
from geopy.distance import vincenty

class drifter_analysis(object):
    
    def __init__(self, data_set=None,time_interval=None):
        
        self.data = data_set
        self.interval = time_interval
        
        self.coordinates = None
        self.velocity = None
        self.analysis()        
        
        # to-be-gemaakte functie voor enigszins relevante data, wellicht. 
        self.time = None
        self.meta_data = None

        
    def analysis(self):
        
        velocities = []
        coordinates = []
        
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
        
            velocities.append((u,v))
            coordinates.append((x0,y0))
            
            if k == len(self.data)//self.interval-1:
                coordinates.append((x1,y1)) 
        
        self.velocity = velocities
        self.coordinates = coordinates
        
