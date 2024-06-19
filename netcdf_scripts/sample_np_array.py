'''
Created on Mar 21, 2024

@author: Karl
'''
import netCDF4 as nc
import numpy as np

def netcdf_to_array_dict(in_ncdf, in_var_list, x_dim, y_dim, dim_dict, preserve_masked_vals=False):
    
    # Create dictionary for ouptut.
    out_dict = {}
    
    # Create dataset.
    dataset = nc.Dataset(in_ncdf, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    dataset.set_auto_mask(preserve_masked_vals)
    
    # Loop through variables in input list.
    for in_var in in_var_list:
    
        # Create variable for variable.
        out_var = dataset[in_var]
        
        # Get dimensions for variable.
        out_var_dims = out_var.dimensions
        
        # Check that input x and y dimensions are in the dataset.
        if set([x_dim, y_dim]).issubset(set(out_var_dims)) == False:
            
            raise Exception("'x_dim' and 'y_dim' are not dimensions in the input netCDF file.")
            
        # Create list to hold modified arrays.
        mod_array_list = [out_var]
        
        # Iterate over dimensions, in order to produce a 2D array.
        for var_idx, var_dim in enumerate(out_var_dims):
            
            # Don't mess with x or y dimensions.
            if var_dim not in [x_dim, y_dim]:
                
                # Loop through dimensions in the input dimension dictionary.
                if var_dim in dim_dict.keys():
                    
                    # Create an index list.
                    var_dim_index_list = list(range(0,dataset.dimensions[var_dim].size))
                    
                    if len(var_dim_index_list) > 1:
                        # Make sure that the value in the dimension dictionary is within the index list.
                        if dim_dict[var_dim] in var_dim_index_list:
                            
                            remove_index_list = var_dim_index_list.remove(dim_dict[var_dim])
                            
                        else:
                            
                            raise Exception("Index '{}' is out of range for dimension '{}'. Check dim_dict input.".format(dim_dict[var_dim], var_dim))

                        # Delete all arrays in the remove_index_list, and squeeze to remove dimension with new length of 1.
                        mod_array = np.squeeze(np.delete(mod_array_list[-1], remove_index_list, axis=var_idx),axis=var_idx)
                        
                        # Add modified array to array list.
                        mod_array_list.append(mod_array)
                    
                else:
                    
                    raise Exception("Dimension '{}' is not in input dim_dict. The dim_dict input must include all non-xy dimensions in the input NetCDF.".format(var_dim))
    
        
        # Add last entry in mod_array_list to output dictionary.
        out_dict[in_var] = np.array(mod_array_list[-1])
        
    # Return dictionary.
    return out_dict


array_dict = netcdf_to_array_dict(r'C:\GIS\Shaw\Waldron_NetCDF\netCDF\TEST_11.nc', ['sub_area', 'easting', 'northing'], 'x', 'y', {'z':5})

print(array_dict)