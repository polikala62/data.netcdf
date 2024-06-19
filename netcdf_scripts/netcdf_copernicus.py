'''
Created on May 14, 2024

@author: Karl
'''

import netCDF4 as nc
import numpy as np
import datetime, ephem

def netcdf_rootgrp(ncdf_path, preserve_masked_vals=False):
    
    # Create dataset.
    rootgrp = nc.Dataset(ncdf_path, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    rootgrp.set_auto_mask(preserve_masked_vals)
    
    # Return rootgrp.
    return rootgrp

def netcdf_var_arrays(rootgrp, lon_var, lat_var, time_var=None):
    
    if time_var == None:
        
        lon_array = np.array(rootgrp.variables[lon_var][:])
        lat_array = np.array(rootgrp.variables[lat_var][:])
 
        return lon_array, lat_array
        
    
    else:
            
        lon_array = np.array(rootgrp.variables[lon_var][:])
        lat_array = np.array(rootgrp.variables[lat_var][:])
        time_array = np.array(rootgrp.variables[time_var][:])
        
        return lon_array, lat_array, time_array


def netcdf_shift_lon(lon_array, in_lon):
    
    # Detect longitude inputs that have been incremented by 360 (i.e. NOAA WCT exports).
    if np.amax(lon_array) > 180: # Longitudes in positive/negative decimal degrees cannot exceed 180.
        return 360 + float(in_lon)
    else:
        return in_lon
    
def netcdf_lonlat_index(lon_array, lat_array, sample_lon, sample_lat):
    
    # If lats and lons are stored as multidimensional arrays, combine them to find indices.
    if len(np.array(lat_array).shape) > 1:
    
        # Subtract select_lat and select_lon from input arrays. 
        lats_diff = np.abs(np.array(lat_array) - float(sample_lat))
        lons_diff = np.abs(np.array(lon_array) - float(sample_lon))
        
        # Find indices corresponding to lat/lon values closest to 'select_lat' and 'select_lon'.
        latlon_diff = np.sqrt((lats_diff **2) + (lons_diff **2))
            
        # Retrieve indices.
        lat_index, lon_index = np.unravel_index(np.argmin(latlon_diff, axis = None), latlon_diff.shape)
    
    # If lats and lons are stored as 1-dimensional arrays (i.e. WCTK exports), search for indices individually.
    else:
        # Find indices corresponding to lat/lon values closest to 'select_lat' and 'select_lon'.
        lat_index_pair = np.unravel_index((np.abs(lat_array-float(sample_lat)).argmin()), lat_array.shape)
        lon_index_pair = np.unravel_index((np.abs(lon_array-float(sample_lon)).argmin()), lon_array.shape)
        
        if len(lat_index_pair) > 1:
            lat_index = int(max(lat_index_pair))
            lon_index = int(max(lon_index_pair))
        else:
            lat_index = int(lat_index_pair[0])
            lon_index = int(lon_index_pair[0])

    return lon_index, lat_index


def netcdf_dimension_index(dim_array, sample_dim_val):
    
    # If lats and lons are stored as multidimensional arrays, combine them to find indices.
    if len(np.array(dim_array).shape) > 1:
    
        # Subtract select_lat and select_lon from input arrays. 
        dim_vals_diff = np.abs(np.array(dim_array) - float(sample_dim_val))
            
        # Retrieve indices.
        dim_index = np.unravel_index(np.argmin(dim_vals_diff, axis = None), dim_array.shape)
    
    # If lats and lons are stored as 1-dimensional arrays (i.e. WCTK exports), search for indices individually.
    else:
        # Find indices corresponding to lat/lon values closest to 'select_lat' and 'select_lon'.
        dim_index = np.unravel_index((np.abs(dim_array-float(sample_dim_val)).argmin()), dim_array.shape)

    return dim_index

def netcdf_datetime_to_days_since_epoch(in_datetime, ref_date, units):
    
    epoch = datetime.datetime.strptime(ref_date,'%Y-%m-%d-%H:%M:%S')
    
    time_diff = (in_datetime - epoch).total_seconds()
    
    if units == "SECONDS":
        
        return int(time_diff)
    
    elif units == "MINUTES":
        
        return int(time_diff / 60)
    
    elif units == "HOURS":
        
        return int(time_diff / 3600)
    
    elif units == "DAYS":
        
        return int(time_diff / 86400)
    
    else:
        
        raise Exception("Parameter 'units' must be 'SECONDS', 'MINUTES', 'DAYS', or 'HOURS'.")
    
#===============================================================================
# MARITIME HORIZONS DATA FUNCTIONS.
#===============================================================================

def read_current_copernicus(in_ncdf, sample_lon, sample_lat, sample_datetime):
    
    u_var = "uo"
    v_var = "vo"
    lon_var = "longitude"
    lat_var = "latitude"
    time_var = "time"
    
    # Use function to convert datetime string to hours since epoch.
    mod_sample_time = netcdf_datetime_to_days_since_epoch(sample_datetime, '1970-01-01-00:00:00', "SECONDS")
    
    # Use functions to get indices.
    rootgrp = netcdf_rootgrp(in_ncdf)
    lonvals, latvals, timevals = netcdf_var_arrays(rootgrp, lon_var, lat_var, time_var)
    mod_sample_lon = netcdf_shift_lon(lonvals, sample_lon)
    
    lon_index = netcdf_dimension_index(lonvals, mod_sample_lon)
    lat_index = netcdf_dimension_index(latvals, sample_lat)
    time_index = netcdf_dimension_index(timevals, mod_sample_time)
    
    # Retrieve value from netcdf.
    sample_u_val = rootgrp.variables[u_var][time_index, lat_index, lon_index]
    sample_v_val = rootgrp.variables[v_var][time_index, lat_index, lon_index]
    
    return [float(i) for i in [sample_u_val, sample_v_val]]
'''
print(read_current_copernicus(r"C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Download\Great_Orme_test_01.nc", 
                              -1.080,
                              50.751,
                              datetime.datetime.strptime('2022-01-01-12:00:00','%Y-%m-%d-%H:%M:%S')))
'''
#------------------------------------------------------------------------------ 

def read_wind_dir_copernicus(in_ncdf, sample_lon, sample_lat, sample_datetime):
    
    dir_var = "wdir10"
    lon_var = "longitude"
    lat_var = "latitude"
    time_var = "time"
    
    # Use function to convert datetime string to hours since epoch.
    mod_sample_time = netcdf_datetime_to_days_since_epoch(sample_datetime, '1970-01-01-00:00:00', "SECONDS")
    
    # Use functions to get indices.
    rootgrp = netcdf_rootgrp(in_ncdf)
    lonvals, latvals, timevals = netcdf_var_arrays(rootgrp, lon_var, lat_var, time_var)
    mod_sample_lon = netcdf_shift_lon(lonvals, sample_lon)
    
    lon_index = netcdf_dimension_index(lonvals, mod_sample_lon)[0]
    
    lat_index = netcdf_dimension_index(latvals, sample_lat)[0]
    time_index = netcdf_dimension_index(timevals, mod_sample_time)
    
    # Retrieve value from netcdf.
    sample_wind_dir_val = rootgrp.variables[dir_var][time_index, lat_index, lon_index]
    
    return float(sample_wind_dir_val)

#------------------------------------------------------------------------------ 

def read_wind_speed_copernicus(in_ncdf, sample_lon, sample_lat, sample_datetime):
    
    speed_var = "si10"
    lon_var = "longitude"
    lat_var = "latitude"
    time_var = "time"
    
    # Use function to convert datetime string to hours since epoch.
    mod_sample_time = netcdf_datetime_to_days_since_epoch(sample_datetime, '1970-01-01-00:00:00', "SECONDS")
    
    # Use functions to get indices.
    rootgrp = netcdf_rootgrp(in_ncdf)
    lonvals, latvals, timevals = netcdf_var_arrays(rootgrp, lon_var, lat_var, time_var)
    mod_sample_lon = netcdf_shift_lon(lonvals, sample_lon)
    
    lon_index = netcdf_dimension_index(lonvals, mod_sample_lon)[0]
    
    lat_index = netcdf_dimension_index(latvals, sample_lat)[0]
    time_index = netcdf_dimension_index(timevals, mod_sample_time)
    
    # Retrieve value from netcdf.
    sample_wind_speed_val = rootgrp.variables[speed_var][time_index, lat_index, lon_index]
    
    return float(sample_wind_speed_val)

#------------------------------------------------------------------------------ 

def read_surface_copernicus(in_ncdf, sample_lon, sample_lat, sample_datetime):
    
    surf_var = "zos"
    lon_var = "longitude"
    lat_var = "latitude"
    time_var = "time"
    
    # Use function to convert datetime string to hours since epoch.
    mod_sample_time = netcdf_datetime_to_days_since_epoch(sample_datetime, '1970-01-01-00:00:00', "SECONDS")
    
    # Use functions to get indices.
    rootgrp = netcdf_rootgrp(in_ncdf)
    lonvals, latvals, timevals = netcdf_var_arrays(rootgrp, lon_var, lat_var, time_var)
    mod_sample_lon = netcdf_shift_lon(lonvals, sample_lon)
    
    lon_index = netcdf_dimension_index(lonvals, mod_sample_lon)
    lat_index = netcdf_dimension_index(latvals, sample_lat)
    time_index = netcdf_dimension_index(timevals, mod_sample_time)
    
    # Retrieve value from netcdf.
    sample_surf_val = rootgrp.variables[surf_var][time_index, lat_index, lon_index]
    
    return [float(i) for i in [sample_surf_val]]

#------------------------------------------------------------------------------ 

def read_vis(in_ncdf, sample_x, sample_y, sample_z):
    
    speed_var = "sub_area"
    x_var = "easting"
    y_var = "northing"
    #z_dim = "z"
    
    # Use functions to get indices.
    rootgrp = netcdf_rootgrp(in_ncdf)
    xvals, yvals = netcdf_var_arrays(rootgrp, x_var, y_var)
    
    x_index = netcdf_dimension_index(xvals, sample_x)[0]
    y_index = netcdf_dimension_index(yvals, sample_y)[0]
    
    #------------------------------------------------------------------------------ 
    # AAAARRRRGGGGHHHHH
    def dict_search(in_dict, search_val):
        
        res_key, res_val = min(in_dict.items(), key=lambda x: abs(search_val - x[1])) #@UnusedVariable
    
        return res_key
    
    z_dict = {1:-8,2:-7,3:-6,4:-5,5:-4,6:-3,7:-2,8:-1,9:0,10:1,11:2,12:3,13:4,14:5,15:6,16:7,17:8}
    
    z_index = dict_search(z_dict, sample_z[0])
    
    #------------------------------------------------------------------------------ 
    
    # Retrieve value from netcdf.
    sample_vis_val = rootgrp.variables[speed_var][z_index, y_index, x_index]
    
    return float(sample_vis_val)

#------------------------------------------------------------------------------ 

#===============================================================================
# FIRST TIDE
#===============================================================================

def twilight_time(in_lon, in_lat, in_datetime, minute_increment, twilight_type, horizon_type="CIVIL"):
    
    # Get datetime for 00:00h for datetime.
    midday_day_string = "{}-12-00-00".format(datetime.datetime.strftime(in_datetime, "%Y-%m-%d"))
    midday_datetime = datetime.datetime.strptime(midday_day_string, "%Y-%m-%d-%H-%M-%S")
    
    # Get horizon value.
    horizon_dict = {"CIVIL":-6, "NAUTICAL":-12, "ASTRONOMICAL":-18}
    horizon = horizon_dict[horizon_type]
    
    # Set sun as astronomical body.
    body = ephem.Sun() #@UndefinedVariableFromImport
    
    # Create ephem observer object.
    obs = ephem.Observer()
    obs.horizon = str(horizon)
    obs.pressure = 0
    obs.date = datetime.datetime.strftime(midday_datetime, "%Y/%m/%d %H:%M:%S")
    obs.lat = str(in_lat)
    obs.lon = str(in_lon)
    
    # Get time for next rising or setting.
    if twilight_type == "SETTING":
        twilight_ephem_time = obs.next_setting(body, start=datetime.datetime.strftime(midday_datetime, "%Y/%m/%d %H:%M:%S"), use_center=True)
    elif twilight_type == "RISING":
        twilight_ephem_time = obs.previous_rising(body, start=datetime.datetime.strftime(midday_datetime, "%Y/%m/%d %H:%M:%S"), use_center=True)
    
    # Convert ephem time to datetime.
    twilight_datetime = datetime.datetime.strptime("-".join(twilight_ephem_time.tuple()),"%Y-%m-%d-%H-%M-%S")
    
    # Round second to nearest time interval.
    twilight_datetime += datetime.timedelta(minutes=(minute_increment * 0.5))
    twilight_datetime -= datetime.timedelta(minutes=twilight_datetime.minute % minute_increment,
                                            seconds=twilight_datetime.second,
                                            microseconds=twilight_datetime.microsecond)
    
    # Return datetime.
    return twilight_datetime

def first_high_tide(in_lon, in_lat, in_datetime, in_minute_inc, in_surface_ncdf, horizon_type="CIVIL"):
    
    # Get first time increment after daylight.
    dawn_time = twilight_time(in_lon, in_lat, in_datetime, in_minute_inc, twilight_type="RISING", horizon_type=horizon_type)
    
    # Create list for check_datetimes.
    check_val_list = []
    
    # Iterate through time increments for a 12-hour period.
    for minute_inc in range(0, 721, in_minute_inc):
        
        # Get datetime.
        check_datetime = dawn_time + datetime.timedelta(minutes=minute_inc)
        
        # Get value from NetCDF.
        check_val_3 = read_surface_copernicus(in_surface_ncdf, in_lon, in_lat, check_datetime)
        
        # Get values for ultimate, penultimate values in list.
        check_val_2 = check_val_list[-1]
        check_val_1 = check_val_list[-2]
        
        # If val 2 is a peak, return datetime.
        if check_val_2 >= check_val_1 and check_val_2 >= check_val_3:
            
            # Return datetime for check_val_2 by subtracting the minute increment from the current iterated value.
            return check_datetime - datetime.timedelta(minutes=minute_inc)
        
        # Add datetime to list.
        check_val_list.append(check_val_3)
        
    # If function completes and no value was returned, return none.
    return None