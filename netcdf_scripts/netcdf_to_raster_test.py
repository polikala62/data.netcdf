'''
Created on 13 Jan 2025

@author: Karl Smith
'''

import arcpy
from netcdf_scripts.netcdf_to_raster import netcdf_to_singleband_raster2

netcdf = ""
x_var = ""
y_var = ""
slice_dict = {}
val_array = []
out_raster_path = ""
out_crs = ""

netcdf_to_singleband_raster2(netcdf, x_var, y_var, slice_dict, val_array, out_raster_path, out_crs, x_dim='', y_dim='')