'''
Created on Mar 24, 2024

@author: Karl
'''
import netCDF4 as nc
'''
def slice_netcdf(ncdf_path, sample_var, x_var, y_var, x_val, y_val, dim_idx_dict, preserve_masked_vals=False):
    
    # Create dataset.
    rootgrp = nc.Dataset(ncdf_path, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    rootgrp.set_auto_mask(preserve_masked_vals)
    
    # Create arrays for x, y, and sample var.
    
    # Process dimensions.
    
    
    
    # Get list of variable names from input and check that sample_var is in the list.
    in_var_dict = rootgrp.variables
    
    # Check that variable exists in var_dict.
    if set([sample_var, x_var, y_var]).issubset(in_var_dict.keys()):
        
        # Create arrays for output.
        sample_array = 
        
    # Raise exception.
    else:
        raise Exception("Variables '{}', '{}', and '{}' are not in file '{}'. Found variables '{}'.".format(sample_var, x_var, y_var, ncdf_path, ", ".join(in_var_dict.keys())))
        
        # Get list of dimensions from input and check them against the dim_dict.
        in_dim_dict = rootgrp.dimensions
        
        # Check that keys in input dictionary exist in the NetCDF dimensions.
        if set(dim_idx_dict.keys()).issubset(in_dim_dict.keys()):
            
            # Create index dictionary, where key is dim name and value is sample index.
            dim_idx_dict = {}
            
            # Add values to index dictionary.
            for in_dim_idx, in_dim in enumerate(in_dim_dict.keys()):
                dim_idx_dict[in_dim] = in_dim_idx
            
            # Create index list, where all indices are 0.
            idx_list = [0 for i in dim_idx_dict.keys()]
            
            # Loop through dimensions, retreiving indices for sample values.
            for sample_dim in dim_dict.keys():
                
                sample_val = dim_dict[sample_dim]
        
        # Raise exception.
        else:
            raise Exception("Dimensions '{}' are not all in file '{}'. Found dimensions '{}'.".format(", ".join(dim_dict.keys()), ncdf_path, ", ".join(in_dim_dict.keys())))
            
    # Raise exception.
    else:
        raise Exception("Variable '{}' is not in file '{}'. Found variables '{}'.".format(sample_var, ncdf_path, ", ".join(in_var_dict.keys())))
    
    # Search arrays for indices, and retrieve values.
'''
def sample_netcdf(ncdf_path, sample_var, var_val_dict, dim_idx_dict, preserve_masked_vals=False): # dim_dict in format {dim_name:sample_val, ...}
    
    # Create dataset.
    rootgrp = nc.Dataset(ncdf_path, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    rootgrp.set_auto_mask(preserve_masked_vals)
    
    # Get list of variable names from input and check that sample_var is in the list.
    in_var_dict = rootgrp.variables
    
    # Check that variable exists in var_dict.
    if sample_var in in_var_dict.keys():
        
        # Get list of dimensions from input and check them against the dim_dict.
        in_dim_dict = rootgrp.dimensions
        
        # Check that keys in input dictionary exist in the NetCDF dimensions.
        if set(dim_idx_dict.keys()).issubset(in_dim_dict.keys()):
            
            # Create index dictionary, where key is dim name and value is sample index.
            dim_idx_dict = {}
            
            # Add values to index dictionary.
            for in_dim_idx, in_dim in enumerate(in_dim_dict.keys()):
                dim_idx_dict[in_dim] = in_dim_idx
            
            # Create index list, where all indices are 0.
            idx_list = [0 for i in dim_idx_dict.keys()]
            
            # Loop through dimensions, retreiving indices for sample values.
            for sample_dim in dim_dict.keys():
                
                sample_val = dim_dict[sample_dim]
        
        # Raise exception.
        else:
            raise Exception("Dimensions '{}' are not all in file '{}'. Found dimensions '{}'.".format(", ".join(dim_dict.keys()), ncdf_path, ", ".join(in_dim_dict.keys())))
            
    # Raise exception.
    else:
        raise Exception("Variable '{}' is not in file '{}'. Found variables '{}'.".format(sample_var, ncdf_path, ", ".join(in_var_dict.keys())))