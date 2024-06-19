'''
Created on Mar 31, 2024

@author: Karl
'''

def add_empty_variable(nc_path, null_val, sample_array, var_name, var_datatype, var_dimensions, var_att_dict):
    
    with nc.Dataset(nc_path, 'w', format="NETCDF4") as dst:
        
        # Create variable.
        dst.createVariable(var_name, var_datatype, var_dimensions)
        
        # Copy variable attributes all at once via dictionary
        dst[var_name].setncatts(var_att_dict)
        
        # Get shape from sample array.
        array_shape = sample_array.shape
        
        # Create empty numpy array.
        array_empty = np.empty(array_shape)
        
        # Fill array with null value.
        array_fill = array_empty.fill(null_val)
        
        # Write array to netcdf.
        dst[var_name][:] = array_fill