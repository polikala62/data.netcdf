'''
Created on Mar 29, 2024

@author: Karl
'''
import arcpy
import netCDF4 as nc
import numpy as np

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
    
pts = r"C:\GIS\CAA_2024\Test_Data\test_coast_poly.shp"
out_crs = arcpy.Describe(pts).spatialReference
out_crs = arcpy.SpatialReference(4326)

netcdf_to_singleband_raster(r'C:\GIS\Maritime_Encounters\OrmeSim\NetCDF\Orme_Surface_2018_01_01.nc',
                            'zos',
                            'longitude',
                            'latitude',
                            'longitude',
                            'latitude',
                            r'C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Raster\copernicus_r_02.tif',
                            out_crs,
                            method="FIRST",
                            flip_v=True)
'''
netcdf_to_singleband_raster(r'C:\GIS\Maritime_Encounters\OrmeSim\Output\full_area_03.nc',
                            'sub_area',
                            'easting',
                            'northing',
                            'x',
                            'y',
                            r'C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Raster\fullarea03_subarea_first.tif',
                            out_crs,
                            method="FIRST",
                            flip_v=True)
'''