import numpy as np
from geopy.distance import vincenty
from datetime import date

class drifter_analysis(object):
    """
    set = True to call certain methods.
      
    """
        
    def __init__(self, data_set=None,time_interval=12, \
                 era_int=False,era_5=False,angle_calc=False):
        
        """      
        :param data_set: list of buoy GPS coordinates. Entries of the form
                        index 0: longitude
                              1: latitude
                              2: year
                              3: month
                              4: day
                              5: hour
                              
        :param time_interval: buoy data is recorded in hourly intervals. Velocity can be 
                         calculated using different time intervals. For example,
                         time_interval = 12 calculates velocity from 12 hour interval changes
                         in location. Possible time_intervals of 12, 24 or 48 hours.
        """
        
        # intialize parameters
        self.data = data_set
        self.interval = int(time_interval)
        
        # calculate coordinates and velocity data using function 'analysis()'
        self.coordinates = None       # list of consecutive GPS coordinates
        self.full_coordinates = None
        self.velocity = None          # list of consecutive velocities (u,v)
        self.speed = None             # absolute velocities
        self.time_set = None          # start time data sampling
        self.start_day = None         # start day data sampling
        self.start_time = None        # full time+day start time (for GRIB use)
        self.days = None              # number of days for which the data runs
        self.analysis()    
        
        # ERA-interim data to be used in angle calculation
        self.u_era_int = None         # u-velocities determined from ERA-Interim
        self.v_era_int = None         # v-velocities determined from ERA-Interim
        if era_int:
            self.era_interim()
            
        # ERA-5 Analysis
        self.u_era_5 = None           # u-velocities determined from ERA-5
        self.v_era_5 = None           # v-velocities determined from ERA-Interim
        if era_5:
            self.era_5_analysis()
        
        self.angles_era_5 = None      # angles between era-5 wind and buoys
        self.angles_era_int = None    # angles between era-interim wind and buoys
        if angle_calc:
            self.angle_calculator(era_int,era_5)
        
    def analysis(self):
        
        """        
        INPUT: buoy GPS data. 
        
        Convert data from GPS-measurements to a list of velocities determined from
        'time_interval' hours of buoy movement. Determined through distance travelled per 
        day/time (forward differencing). 
        
        Important: data should be trimmed to full days, 
        i.e. the data starts at 00:00 of day 1 and ends at 23 hours of the final day. 
            
        OUTPUT: Returns (sets self.parametes) list with consecutive velocities between period (hrs) of set 
        interval, as well as meta-deta concerning the buoy movement (coordinates,start-time, duration, absolute velocity).   
           
        """
        
        # Determine the number of days in the data set
        date_1 = date(int(self.data[0][2]),int(self.data[0][3]),int(self.data[0][4]))
        date_2 = date(int(self.data[-1][2]),int(self.data[-1][3]),int(self.data[-1][4]))
        data_diff = date_2 - date_1
        
        self.start_day = int(self.data[0][4])   
        self.days = int(data_diff.days) + 1                               # + 1 because date package calculates 'in between' days  
        
        velocities = np.zeros((int(self.days*24/self.interval)-1,2))      # N-1 velocity data points using central or forward approximation
        speed = np.zeros(int(self.days*24/self.interval)-1)               # absolute velocities
        coordinates = np.zeros((int(self.days*24/self.interval),2))       # coordinates (lon, lat) at time intervals
        time_msr = np.zeros((24,0)).tolist()                              # store indices of measurements at fixed times of day

        for measurement in self.data:
            time = measurement[-1]
            time_msr[int(time)].append(measurement)
            
        for element in range(24):
            if len(time_msr[element]) == self.days:
                if len(time_msr[element + 6]) == self.days:
                    if len(time_msr[element + 12]) == self.days:
                        if len(time_msr[element + 18]) == self.days:
                            self.time_set = element
                            break
                        
        # error message if the data is not suitable for 6 hour intervals
        if self.time_set == None:
            print('No continuous sequence of measurements at chosen hour and interval.')
        
        # Create lists of consecutive latitudes and longitudes of measurements   
        lons = np.zeros(int(self.days*4))
        lats = np.zeros(int(self.days*4))
       
        # populate lons and lats with 6-hourly interval data
        for k in range(self.days):
            lons[k*4] = time_msr[self.time_set][k][0] 
            lons[k*4+1] = time_msr[self.time_set+6][k][0]
            lons[k*4+2] = time_msr[self.time_set+12][k][0]
            lons[k*4+3] = time_msr[self.time_set+18][k][0]
            
            lats[k*4] = time_msr[self.time_set][k][1] 
            lats[k*4+1] = time_msr[self.time_set+6][k][1]
            lats[k*4+2] = time_msr[self.time_set+12][k][1]
            lats[k*4+3] = time_msr[self.time_set+18][k][1]
            
        full_coordinates = np.zeros((len(lons),2))
        
        for m in range(len(lons)):
            full_coordinates[m] = [lons[m],lats[m]]
            
        self.full_coordinates = full_coordinates
            
        # trim the data to match the interval (input)
        lons = lons[::int(self.interval/6)]
        lats = lats[::int(self.interval/6)]
        
        # calculate velocities (average per day: first entry is velocity on day 1) between self.interval points
        for k in range(1,len(lons)):          
            x0 = lons[k-1]
            x1 = lons[k]
            y0 = lats[k-1]
            y1 = lats[k]
       
            # vincenty is a method to determine distance between points on earth. Longitudonal distance for u, latitudonal for v.
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
            if k == (len(lons)-1):
                coordinates[k] = [x1,y1] 
              
        self.velocity = velocities
        self.coordinates = coordinates
        self.speed = speed
        
        # of the form 201312010 for use in pygrib 'keys' function call
        self.start_time = str(date_1.year)+str(date_1.month)+str(date_1.day)+1*str(0)
        
    def era_interim(self):
               
        """
        lorem ipsum
        
        """
        import pygrib 
        
        # fixed lat/lon grid of interest for drifters
        lat1 = -70      # South
        lat2 = -40      # North
        lon1 = 140      # West
        lon2 = 230      # East
               
        era_int = pygrib.open('/home/wim/Documents/MAIO Drifters/era_interim/ERA_Interim.grib')
   
        # trim the data to the right start time (excluding off set due to potentially missing data)
        counter = 0        
        
        for grb in era_int:
            if str(grb.dataDate)+str(grb.dataTime) == self.start_time:
                break
            else:
                counter += 1
                
        era_int_u = era_int.select(name='10 metre U wind component')
        era_int_v = era_int.select(name='10 metre V wind component')
        era_int.close()
        print('ERA-interim counter is set to '+str(counter)+'. This means that the first '+str(int(counter/2/4))+' days of the data set are skipped.')
        
        # the actual trim using the off set 'counter'. Divide counter by two due to u and v presence. 
        era_int_u = era_int_u[counter//2:]
        era_int_v = era_int_v[counter//2:]        
        
        # create a full list of all the era-interim data at the points where velocities are determined.
        # For each velocity measurement, era requires 2, 4 and 8 data-points for 12, 24 and 48 time interval, respectively.
        full_era_u = np.zeros(len(self.full_coordinates))
        full_era_v = np.zeros(len(self.full_coordinates))

        for l in range(len(self.full_coordinates)):
            
            # these are the coordinates we are interested in linking to the velocity field
            lon_q = self.full_coordinates[l][0]
            lat_q = self.full_coordinates[l][1]
                                       
            data_u,lat,lon = era_int_u[l].data(lat1=lat1,lat2=lat2,lon1=lon1,lon2=lon2)
            data_v,lat,lon = era_int_v[l].data(lat1=lat1,lat2=lat2,lon1=lon1,lon2=lon2)
            lat_count = 0
            lon_count = 0
                
            # lat runs from ~90 to -90
            for j in lat[:,0]:
                if j < lat_q:
                    break
                else:
                    lat_count += 1
            # lon runs from 0 to ~360
            for j in lon[0,:]:
                if j > lon_q:
                    break
                else:
                    lon_count += 1
                
            # check which neighbouring lon/lat is closest to buoy location
            if np.abs(lat[lat_count+1,0] - lat_q) < np.abs(lat[lat_count,0] - lat_q):
                lat_count += 1
            if np.abs(lon[0,lon_count+1] - lat_q) < np.abs(lon[0,lon_count] - lat_q):
                lon_count += 1
                
            full_era_u[l] = data_u[lat_count][lon_count]
            full_era_v[l] = data_v[lat_count][lon_count]
        
        era_interim_u = np.zeros(len(self.velocity))
        era_interim_v = np.zeros(len(self.velocity))
        
        n = self.interval//6
        temp_u = [ full_era_u[i:i+n] for i in range(0,len(full_era_u),n)]
        temp_v = [ full_era_v[i:i+n] for i in range(0,len(full_era_v),n)]
        
        for k in range(len(self.velocity)):
            era_interim_u[k] = np.mean(temp_u[k])
            era_interim_v[k] = np.mean(temp_v[k])

        self.u_era_int = era_interim_u
        self.v_era_int = era_interim_v

    def era_5_analysis(self):
            
        """       
        This method relies on a numpy array of the era_5 data as well as longitude and latitude 'mapping' arrays.
        
        INPUT: numpy array of dimension (latitude, longitude, measurements (t)) for wind field, and latitude/longitude for 
        the 'mapping', i.e. the lat/long grid identifies the location of the wind field data points.
        
        OUTPUT: wind velocities (u,v) for requested time interval at GPS locations as determined by era-5 reanalysis data. 
        
        """
                           
        lat = np.load('/home/wim/Documents/MAIO Drifters/era_5/lat_era_5_data-70_-40_140_230.npy')
        lon = np.load('/home/wim/Documents/MAIO Drifters/era_5/lon_era_5_data-70_-40_140_230.npy')
                      
        data_u = np.load('/home/wim/Documents/MAIO Drifters/era_5/np_era_5_u.npy')
        data_v = np.load('/home/wim/Documents/MAIO Drifters/era_5/np_era_5_v.npy')

        # trim data to match start time buoy data
        day_off_set = self.start_day - 1   # The -1 has been checked and double checked asnd is indeed correct/required
        
        data_u = data_u[:,:,int(24*day_off_set+self.time_set)::6]
        data_v = data_v[:,:,int(24*day_off_set+self.time_set)::6]

        era_5_u_full = np.zeros(len(self.full_coordinates))
        era_5_v_full = np.zeros(len(self.full_coordinates))
        
        # sample every 6th (u,v)
        for k in range(len(self.full_coordinates)):
    
            lon_q = self.full_coordinates[k][0]
            lat_q = self.full_coordinates[k][1]

            lat_count = 0
            lon_count = 0
                
            # lat runs from ~90 to -90
            for j in lat[:,0]:
                if j < lat_q:
                    break
                else:
                    lat_count += 1
            # lon runs from 0 to ~360
            for j in lon[0,:]:
                if j > lon_q:
                    break
                else:
                    lon_count += 1
             
            # check which neighbouring lon/lat is closest to buoy location
            if np.abs(lat[lat_count+1,0] - lat_q) < np.abs(lat[lat_count,0] - lat_q):
                lat_count += 1
            if np.abs(lon[0,lon_count+1] - lat_q) < np.abs(lon[0,lon_count] - lat_q):
                lon_count += 1
            
            era_5_u_full[k] = data_u[lat_count,lon_count,k]
            era_5_v_full[k] = data_v[lat_count,lon_count,k]
            
        era_5_u = np.zeros(len(self.velocity))
        era_5_v = np.zeros(len(self.velocity))
        
        n = self.interval//6
        temp_u = [ era_5_u_full[i:i+n] for i in range(0,len(era_5_u_full),n)]
        temp_v = [ era_5_v_full[i:i+n] for i in range(0,len(era_5_v_full),n)]
       
        for k in range(len(self.velocity)):
            era_5_u[k] = np.mean(temp_u[k])
            era_5_v[k] = np.mean(temp_v[k])

        self.u_era_5 = era_5_u
        self.v_era_5 = era_5_v
        print(self.time_set)
        
    def angle_calculator(self,era_int,era_5):
            
        """
        This method calculates the angle between the buoy velocity and either era-5 or era-interim data.
        
        INPUT: buoy, era-5 and era-interim velocity fields as determined by the analysis methods present in this class.            
        OUTPUT: List of consecutive angles between buoy and era-5/era-interim wind field.
            
        """
            
        if era_5:
            angles_era_5 = np.zeros(len(self.velocity))
            
            for k in range(len(self.velocity)):
                v0 = [self.u_era_5[k],self.v_era_5[k]]
                v1 = self.velocity[k]
                angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
                # Positive angle --> buoy movement to the right of the wind
                angles_era_5[k] = np.degrees(angle)
                    
            self.angles_era_5 = angles_era_5
                    
        if era_int:
            angles_era_int = np.zeros(len(self.velocity))
                
            for k in range(len(self.velocity)):
                v0 = [self.u_era_int[k],self.v_era_int[k]]
                v1 = self.velocity[k]
                angle = -np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
                # Positive angle --> buoy movement to the right of the wind
                angles_era_int[k] = np.degrees(angle)
               
            self.angles_era_int = angles_era_int