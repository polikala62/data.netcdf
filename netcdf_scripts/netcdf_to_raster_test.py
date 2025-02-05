'''
Created on 13 Jan 2025

@author: Karl Smith
'''

import arcpy
import netCDF4 as nc
import numpy as np
from netcdf_scripts.netcdf_to_raster import netcdf_to_singleband_raster2
from netcdf_scripts import numpy_array

netcdf = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\malta_test_34.nc"
x_var = "easting"
y_var = "northing"
slice_dict = {"x":[0,529],
              "y":[0,237],
              "z":[0,1],
              "d":[42,43]}


val_slice = numpy_array.slice_netcdf(netcdf, "sub_area", slice_dict, preserve_masked_vals=True)
val_array = numpy_array.flatten_slice_by_dict(netcdf, "sub_area", val_slice, {'z':'MEAN', 'd':'MEAN'})

out_raster_path = r"C:\GIS\ArcPro_Projects\Visibility_Testing\Rasters\malta_test_34.tif"
#out_crs = arcpy.SpatialReference(4326)
out_crs = arcpy.Describe(r'C:\GIS\ArcPro_Projects\PM_Med_Survey\Features\med_vis\inputs\vis_buffer_04.shp').spatialReference
x_dim = "x"
y_dim = "y"

netcdf_to_singleband_raster2(netcdf, x_var, y_var, slice_dict, val_array, out_raster_path, out_crs, x_dim, y_dim)