'''
Created on Mar 29, 2024

@author: Karl
'''
import arcpy
import netCDF4 as nc
import numpy as np

from netcdf_scripts import numpy_array

def netcdf_to_singleband_raster2(netcdf, x_var, y_var, slice_dict, val_array, out_raster_path, out_crs, x_dim='', y_dim=''):
    
    # Set arcpy environments.
    arcpy.env.outputCoordinateSystem = out_crs
    arcpy.env.overwriteOutput = True
    
    # Set x and y dimensions to the same as the variables (they usually are), unless values have been supplied.
    if x_dim == '':
        x_dim = x_var
    if y_dim == '':
        y_dim = y_var
    
    # Generate arrays for x and y variables.
    x_array = numpy_array.slice_netcdf(netcdf, x_var, slice_dict)
    y_array = numpy_array.slice_netcdf(netcdf, y_var, slice_dict)
    
    # Flatten x and y arrays.
    flat_x_array = numpy_array.flatten_slice(netcdf, x_var, x_array, "AVERAGE", [x_dim])
    flat_y_array = numpy_array.flatten_slice(netcdf, y_var, y_array, "AVERAGE", [y_dim])
    
    # Calculate average resolution of raster.
    x_res = float(abs(np.average(np.diff(flat_x_array))))
    y_res = float(abs(np.average(np.diff(flat_y_array))))
    
    
    # KLUDGE
    y_res=x_res
    
    # Check that values in x and y arrays are regular. If they aren't, print a warning (and ignore it).
    if (np.diff(flat_x_array)==np.diff(flat_x_array)[0][0]).all() == False:
        print("WARNING: netcdf_to_singleband_raster: values in x variable are not evenly spaced. Using an average spacing of {}.".format(x_res))
    if (np.diff(flat_y_array)==np.diff(flat_y_array)[0][0]).all() == False:
        print("WARNING: netcdf_to_singleband_raster: values in y variable are not evenly spaced. Using an average spacing of {}.".format(y_res))
    
    # Caluclate extent.
    x_min = float(np.nanmin(flat_x_array))
    x_max = float(np.nanmax(flat_x_array))
    y_min = float(np.nanmin(flat_y_array))
    y_max = float(np.nanmax(flat_y_array))
    
    #arcpy.env.extent = "{} {} {} {}".format(x_min, y_min, x_max, y_max)
    
    # Calculate lower-left corner.
    sw_corner = arcpy.Point(x_min, y_min)
    
    # Create list to hold modified value arrays.
    mod_val_array_list = [val_array]
    # Check orientation of flattened arrays. If array (0,0) is upper-right corner, then values should be equivalent to min(x) and max(y).
    if x_min != flat_x_array[0][0]:
        print("WARNING: netcdf_to_singleband_raster: flipping along vertical axis...")
        mod_val_array_list.append(np.flipLR(mod_val_array_list[-1]))
    print(y_max, flat_y_array[0][0])
    if y_max != flat_y_array[0][0]:
        print("WARNING: netcdf_to_singleband_raster: flipping along horizontal axis...")
        mod_val_array_list.append(np.flipud(mod_val_array_list[-1]))
        
    # Replace mask values with arbitrary value.
    fill_val = 999999
    #mod_val_array_list.append(mod_val_array_list[-1].filled(fill_val))
    mod_val_array_list.append(np.ma.filled(mod_val_array_list[-1], fill_val))
    
    # Replace really high values with None (TODO: CHANGE THIS).
    mod_val_array_list.append(np.where(mod_val_array_list[-1] < fill_val, mod_val_array_list[-1], fill_val))
    
    # Create and save raster.
    out_ras = arcpy.NumPyArrayToRaster(mod_val_array_list[-1], lower_left_corner=sw_corner, x_cell_size=x_res, y_cell_size=y_res, value_to_nodata=fill_val)
    out_ras.save(out_raster_path)
    

def netcdf_to_singleband_raster(in_ncdf, in_var, x_var, y_var, x_dim, y_dim, out_raster, out_crs, 
                                method="AVERAGE", flip_v=False, flip_h=False, preserve_masked_vals=False):
    
    # Set arcpy environments.
    arcpy.env.outputCoordinateSystem = out_crs
    arcpy.env.overwriteOutput = True
    
    # Create dataset.
    dataset = nc.Dataset(in_ncdf, "r", format="NETCDF4")
    
    # Get rid of masked values (!!!)
    dataset.set_auto_mask(preserve_masked_vals)
    
    # Create variable for variable.
    out_var = dataset[in_var]
    
    # Get dimensions for variable.
    out_var_dims = out_var.dimensions
    
    # Check that x and y dimensions are in the dataset.
    if x_dim in out_var_dims and y_dim in out_var_dims:
        
        # Create list to hold modified arrays.
        mod_array_list = [out_var]
        
        # Iterate over dimensions, in order to produce a 2D array.
        for var_idx, var_dim in enumerate(out_var_dims):
            
            # Don't mess with x or y dimensions.
            if var_dim not in [x_dim, y_dim]:
            
                if method == "AVERAGE":
                    
                    # Average the last entry in mod_array_list by dimension.
                    mod_array = np.average(mod_array_list[-1],axis=var_idx)
                    
                    # Add averaged array to list.
                    mod_array_list.append(mod_array)
                    
                elif method == "FIRST":
                    
                    # Create an index list that starts at 1 and increases to the length of the dimension.
                    remove_indices = range(1,dataset.dimensions[var_dim].size)
                    
                    if len(remove_indices) > 1:
                        
                        # Delete all arrays after the first, and squeeze to remove dimension with new length of 1.
                        mod_array = np.squeeze(np.delete(mod_array_list[-1], remove_indices, axis=var_idx),axis=var_idx)
                    
                        # Add modified array to array list.
                        mod_array_list.append(mod_array)
                    
                elif method == "LAST":
                    
                    # Create an index list that starts at 1 and increases to the length of the dimension.
                    remove_indices = range(0,(dataset.dimensions[var_dim].size-1))
                    
                    # Delete all arrays after the first, and squeeze to remove dimension with new length of 1.
                    mod_array = np.squeeze(np.delete(mod_array_list[-1], remove_indices, axis=var_idx),axis=var_idx)
                    
                    # Add modified array to array list.
                    mod_array_list.append(mod_array)
                    
                else:
                    
                    raise Exception("Method '{}' is not supported.".format(method))
    
    else:
        
        raise Exception("'x_dim' and 'y_dim' are not dimensions in the input netCDF file.")
    
    # FIND RESOLUTION
    x_min, y_min = [np.min(dataset[i]) for i in (x_var, y_var)]
    
    x_max, y_max = [np.max(dataset[i]) for i in (x_var, y_var)]
    
    x_range, y_range = [dataset.dimensions[i].size for i in (x_dim, y_dim)]
    
    x_res = (x_max - x_min) / x_range
    y_res = (y_max - y_min) / y_range
    
    # FIND CORNER POINT
    sw_corner = arcpy.Point(float(x_min), float(y_min))
    
    # CREATE RASTER
    
    if flip_v:
        
        mod_array_list.append(np.flip(mod_array_list[-1],0))
        
    if flip_h:
        
        mod_array_list.append(np.flip(mod_array_list[-1],1))
        
    out_array = mod_array_list[-1]
    
    out_ras = arcpy.NumPyArrayToRaster(out_array, lower_left_corner=sw_corner, x_cell_size=x_res, y_cell_size=y_res, value_to_nodata=None)
    
    out_ras.save(out_raster)
    
#pts = r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp"
#out_crs = arcpy.Describe(pts).spatialReference
#out_crs = arcpy.SpatialReference(4326)
#wgs = arcpy.SpatialReference("WGS 1984")
#in_netcdf = r"C:\GIS\NetCDF_Manipulation\Data\test_echam_spectral.nc"
#slice_dict = {'time':[0,8], 'mlev':[0,47], 'lat':[0,96], 'lon':[0,192]}
#pr_slice = numpy_array.slice_netcdf(in_netcdf, 'xl', slice_dict)
#slice_avg = numpy_array.flatten_slice_by_dict(in_netcdf, 'xl', pr_slice, {'time':"MEAN", 'mlev':"MEAN"})
#netcdf_to_singleband_raster2(in_netcdf, 'lon', 'lat', slice_dict, slice_avg, r'C:\GIS\NetCDF_Manipulation\Output\test_echam_01.tif', wgs)
