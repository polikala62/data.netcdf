'''
Created on 4 Nov 2024

@author: Karl Smith
'''

import arcpy
import netCDF4 as nc
import numpy as np

# Flattens array around input axis.
def filter_netcdf_array(netcdf, variable, in_slice, method, exclude_dims):
    
    #TODO: Add any kind of validation.
    
    # Create dataset.
    dataset = nc.Dataset(netcdf, "r", format="NETCDF4")
    
    # Get netCDF variable as python variable (I KNOW)
    var = dataset.variables[variable]
    
    # Get list of dimensions from the variable.
    dim_list = list(var.dimensions)