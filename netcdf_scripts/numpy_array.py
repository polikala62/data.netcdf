'''
Created on 21 Oct 2024

@author: Karl Smith
'''

import netCDF4 as nc
import numpy as np

#------------------------------------------------------------------------------ 

# Prints dimensions and shape for variable, to help generate slice_dict (see below).
def variable_dimensions_list(netcdf, variable):
    
    # Create dataset.
    rootgrp = nc.Dataset(netcdf, "r", format="NETCDF4")
    var = rootgrp.variables[variable]
    
    return list(zip(var.dimensions, var.shape))
    
#------------------------------------------------------------------------------ 

# Indexes ranges in a dictionary a numpy array from a netCDF file.
def slice_netcdf(netcdf, variable, slice_dict, preserve_masked_vals=True): #slice_dict should be in format {'dimension name':[start index, stop index]}
    
    # Create dataset.
    dataset = nc.Dataset(netcdf, "r", format="NETCDF4")
    
    # Get rid of masked values (or don't, unmasked values will screw up stats)
    dataset.set_auto_mask(preserve_masked_vals)
    
    # Get netCDF variable as python variable (I KNOW)
    var = dataset.variables[variable]
    
    # Get list of dimensions from the variable.
    dim_list = var.dimensions
    
    # Check dimensions in input dictionary against dimensions attached to variable.
    
    # If elements in input dictionary are not in netcdf, print warning and ignore.
    if len(list(set(slice_dict.keys()) - set(dim_list))) > 0:
        warning_str = "'{}'".format("', '".join(list(set(slice_dict.keys()) - set(dim_list))))
        print("WARNING: Dimension(s) {} could not be found in input netCDF file. Ignoring {}.".format(warning_str, warning_str))
        
    # If elements are in netcdf are not in dictionary, print warning and abort.
    if (len(list(set(dim_list) - set(slice_dict.keys())))) > 0:
        warning_str = "'{}'".format("', '".join(list(set(dim_list) - set(slice_dict.keys()))))
        raise Exception("Dimension(s) {} are in the netCDF file but not in slice_dict. Could not generate slice.".format(warning_str))
    
    # Check that all input dictionary entries are list objects with length of 2.
    try:
        if False in [len(slice_dict[i]) == 2 for i in slice_dict.keys()]:
            raise Exception("Values in 'slice_dict' must be list objects with two items, i.e. [bottom_bound, top_bound].")
    except:
        raise Exception("Values in 'slice_dict' must be list objects with two items, i.e. [bottom_bound, top_bound].")
    
    # VALIDATION TO ADD: Values can't be out of bounds, can't be negative, lower bound can't be greater than upper bound.
    
    # Create list to hold slice objects.
    slice_obj_list = [slice(slice_dict[i][0], slice_dict[i][1]) for i in list(dim_list)]
    
    # Return sliced variable as array.
    return var[slice_obj_list]

#------------------------------------------------------------------------------ 

# Flattens array around input axis.
def flatten_slice(netcdf, variable, in_slice, method, exclude_dims):
    
    #TODO: Add any kind of validation.
    
    # Create dataset.
    dataset = nc.Dataset(netcdf, "r", format="NETCDF4")
    
    # Get netCDF variable as python variable (I KNOW)
    var = dataset.variables[variable]
    
    # Get list of dimensions from the variable.
    dim_list = list(var.dimensions)
    
    # Create list of arrays.
    pr_array_list = [in_slice]
    
    # Remove excluded dimensions from list.
    mod_dim_list = [i for i in dim_list if i not in exclude_dims]
    
    # Loop through items in flatten_dict.
    for dim_name in mod_dim_list:
    
        # Find index for input dimension.
        dim_idx = list(dim_list).index(dim_name)
        
        # Add flattened array to processing list.
        if method == "MAXIMUM":
            
            # Return slice with maximum values.
            pr_array_list.append(np.nanmax(pr_array_list[-1], axis=dim_idx))
        
        elif method == "MEAN":
        
            # Return slice averaged around input dimension.
            pr_array_list.append(np.nanmean(pr_array_list[-1], axis=dim_idx))
        
        elif method == "MINIMUM":
            
            # Return slice with maximum values.
            pr_array_list.append(np.nanmin(pr_array_list[-1], axis=dim_idx))
        
        elif method == "VARIANCE":
            
            # Return slice with variance values.
            pr_array_list.append(np.nanvar(pr_array_list[-1], axis=dim_idx))
        
        elif method == "STANDARD_DEVIATION":
            
            # Return slice with variance values.
            pr_array_list.append(np.nanstd(pr_array_list[-1], axis=dim_idx))
            
    return pr_array_list[-1]

# Flattens array around input axis using a dictionary to select statistic.
def flatten_slice_by_dict(netcdf, variable, in_slice, flatten_dict):
    
    #TODO: Add any kind of validation.
    
    # Create dataset.
    dataset = nc.Dataset(netcdf, "r", format="NETCDF4")
    
    # Get netCDF variable as python variable (I KNOW)
    var = dataset.variables[variable]
    
    # Get list of dimensions from the variable.
    dim_list = list(var.dimensions)
    
    # Create list of arrays.
    pr_array_list = [in_slice]
    
    # Loop through items in flatten_dict.
    for dim_name in flatten_dict.keys():
        
        method = flatten_dict[dim_name]
    
        # Find index for input dimension.
        dim_idx = list(dim_list).index(dim_name)
    
        # Add flattened array to processing list.
        if method == "MAXIMUM":
            
            # Return slice with maximum values.
            pr_array_list.append(np.nanmax(pr_array_list[-1], axis=dim_idx))
        
        elif method == "MEAN":
        
            # Return slice averaged around input dimension.
            pr_array_list.append(np.nanmean(pr_array_list[-1], axis=dim_idx))
        
        elif method == "MINIMUM":
            
            # Return slice with maximum values.
            pr_array_list.append(np.nanmin(pr_array_list[-1], axis=dim_idx))
        
        elif method == "VARIANCE":
            
            # Return slice with variance values.
            pr_array_list.append(np.nanvar(pr_array_list[-1], axis=dim_idx))
            
        # Update dim_list.
        dim_list = [i for i in dim_list if i != dim_name]
        print(dim_list)
        
    return pr_array_list[-1]
    
#------------------------------------------------------------------------------ 
pass

#in_netcdf = r"C:\GIS\NetCDF_Manipulation\Data\test_echam_spectral.nc"
#pr_slice = slice_netcdf(in_netcdf, 'xl', {'time':[0,8], 'mlev':[0,47], 'lat':[0,96], 'lon':[0,192]})
#print(pr_slice.shape)
#slice_avg = flatten_slice(in_netcdf, 'xl', pr_slice, {'time':"MEAN", 'mlev':"MEAN"})
#print(slice_avg.shape)

#print(variable_dimensions_list(in_netcdf, 'xl'))
