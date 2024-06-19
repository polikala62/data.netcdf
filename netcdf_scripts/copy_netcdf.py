'''
Created on Mar 31, 2024

@author: Karl
'''

def copy_netcdf(in_ncdf_path, out_ncdf_path, var_list): # Field dict should be in format [field_name]:[upper bound, lower bound]
    
    with nc.Dataset(in_ncdf_path, 'r', format="NETCDF4") as src, nc.Dataset(out_ncdf_path, 'w', format="NETCDF4") as dst:
        # copy global attributes all at once via dictionary
        dst.setncatts(src.__dict__)
        
        # copy dimensions
        for name, dimension in src.dimensions.items():
            #print(name)
            #print(dimension)
            #print()
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if name in var_list:
                #print(name)
                #print(variable)
                print(variable.datatype)
                
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)
                #print(src[name].__dict__)
                dst[name][:] = src[name][:]
                