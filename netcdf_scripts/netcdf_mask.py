
import arcpy
import numpy as np
from netcdf import netcdf_copernicus as netcdf
import datetime

def netcdf_mask(in_extent_fc, xy_spacing, path_dict, mask_dict, in_datetime, out_raster, in_coast_polyline=None, coast_dist_limit=None):
    
    def coordinate_distance_to_geometry(coord_x, coord_y, geometry, crs):
        
        coord_geometry = arcpy.PointGeometry(arcpy.Point(coord_x, coord_y), crs)
        
        distance_list = []
        
        # Loop through features in input geometry.
        for feature in geometry:
            
            distance_list.append(coord_geometry.distanceTo(feature))
        
        return min(distance_list)
    
    # Find resolution for output array (taken from vis.radial_viewshed).
    
    def crs_to_WGS84(xy, in_crs): #xy must be tuple: [x, y] in metres.
        
        # Create in-memory feature class.
        temp_fc = r'in_memory\ETRS2WGS'
        
        if arcpy.Exists(temp_fc) == True:
            arcpy.Delete_management(temp_fc)
    
        # Set variables for spatial references.
        WGS84 = arcpy.SpatialReference(4326)
        
        # Create numpy array from input.
        xy_tuple = (xy[0], xy[1])
        array = np.array([(1, xy_tuple)], np.dtype([('idfield',np.int32),('XY', '<f8', 2)]))
    
        # Convert numpy array to feature class, using WGS84.
        arcpy.da.NumPyArrayToFeatureClass(array, temp_fc, ['XY'], in_crs) #@UndefinedVariableFromImport
        
        # Convert feature class to numpy array, using ETRS89.
        res = arcpy.da.FeatureClassToNumPyArray(temp_fc, 'SHAPE@XY', spatial_reference=WGS84) #@UndefinedVariableFromImport
        
        # Return values from array.
        return [res[0][0][0], res[0][0][1]]
    
    #===========================================================================
    # IMPORT LINE FEATURES TO GEOMETRY.
    #===========================================================================
    
    #@TODO.
    if in_coast_polyline != None:
        coastline_features = arcpy.CopyFeatures_management(in_coast_polyline, arcpy.Geometry())
    
    #===========================================================================
    # GET CRS FROM INPUT FEATURES.
    #===========================================================================
    
    extent_fc_crs = arcpy.Describe(in_extent_fc).spatialReference
    
    #===========================================================================
    # GET X, Y RANGES FROM INPUT MASK SHAPE.
    #===========================================================================
    
    # Create an extent object for the land mask.
    land_fc_extent = arcpy.Describe(in_extent_fc).extent
    
    # Round west and south up to the nearest spacing value.
    x_min = land_fc_extent.XMin + (xy_spacing - (land_fc_extent.XMin % xy_spacing))
    y_min = land_fc_extent.YMin + (xy_spacing - (land_fc_extent.YMin % xy_spacing))
    
    # Round east and north down to the nearest spacing value.
    x_max = land_fc_extent.XMax - (land_fc_extent.XMax % xy_spacing)
    y_max = land_fc_extent.YMax - (land_fc_extent.YMax % xy_spacing)
    
    # Calculate ranges.
    #x_range = [land_fc_west, land_fc_east, xy_spacing]
    #y_range = [land_fc_south, land_fc_north, xy_spacing]
    
    
    x_range = np.arange(x_min, (x_max + xy_spacing), xy_spacing)
    y_range = np.arange(y_min, (y_max + xy_spacing), xy_spacing)
    
    # Delete describe object.
    del land_fc_extent
    
    # Create array for mask intersection.
    out_array_shape = (len(y_range), len(x_range))
    iter_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    out_array = np.zeros(out_array_shape, dtype=float, order='C')
    
    # Loop through xy coordinates.
        
    # Iterate through indices.
    for index, x in np.ndenumerate(iter_array): #@UnusedVariable
        
        # Get indices from iterator.
        y_idx, x_idx = index
        
        # Get x, y values for iterated data point.
        x_val = float(x_range[x_idx])
        y_val = float(y_range[y_idx])
        
        # Check point distance to coast polyine.
        xy_coast_dist = coordinate_distance_to_geometry(x_val, y_val, coastline_features, extent_fc_crs)
        
        if in_coast_polyline != None and coast_dist_limit != None:
            
            if xy_coast_dist <= coast_dist_limit:
            
                process_pt = True
            
            else:
                
                process_pt = False
            
        else:
            
            process_pt = True
            
        
        if process_pt:
            
            # Convert coordinates to WGS84.
            #iter_lon, iter_lat = crs_to_WGS84([x_val, y_val], extent_fc_crs)
            iter_lon, iter_lat = x_val, y_val
            
            # Collect sea level from NetCDF.
            real_elevation = netcdf.read_surface_copernicus(path_dict['surface'], iter_lon, iter_lat, in_datetime)
            
            # Create vector for real current direction at point.
            water_velocity_x, water_velocity_y = netcdf.read_current_copernicus(path_dict['water'], iter_lon, iter_lat, in_datetime)
            
            # Create vector for real wind direction at point.
            wind_dir = netcdf.read_wind_dir_copernicus(path_dict['wind_dir'], iter_lon, iter_lat, in_datetime)
            wind_speed = netcdf.read_wind_speed_copernicus(path_dict['wind_speed'], iter_lon, iter_lat, in_datetime)
            
            if mask_dict['surface'] == real_elevation or mask_dict['water_x'] == water_velocity_x or mask_dict['water_y'] == water_velocity_y or mask_dict['wind_dir'] == wind_dir or mask_dict['wind_speed'] == wind_speed:
                
                out_array[y_idx, x_idx] = 0
                
            else:
                
                out_array[y_idx, x_idx] = 1
            
        else:
            
            out_array[y_idx, x_idx] = 2
        
        # Create raster.
        sw_corner = arcpy.Point(float(y_min), float(x_min))
        flip_array = np.flip(out_array,0)
        out_ras = arcpy.NumPyArrayToRaster(flip_array, lower_left_corner=sw_corner, x_cell_size=xy_spacing, y_cell_size=xy_spacing, value_to_nodata=None)
        out_ras.save(out_raster)
    
    print("Script finished!")

#------------------------------------------------------------------------------ 
 
netcdf_mask(r"C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Features\DEM_footprint_WGS84.shp", 
            50000, 
            {'surface':"C:\\GIS\\Maritime_Encounters\\OrmeSim\\NetCDF\\Orme_Surface_2018_01_01.nc",'water':"C:\\GIS\\Maritime_Encounters\\OrmeSim\\NetCDF\\Orme_Water_2018_01_01.nc",'wind_dir':"C:\\GIS\\Maritime_Encounters\\Copernicus_Testing\\Download_Test\\Download\\adaptor.mars.external-1715692160.536764-3213-16-1487bebf-9495-410f-8dba-2e88baf99de9.nc",'wind_speed':"C:\\GIS\\Maritime_Encounters\\Copernicus_Testing\\Download_Test\\Download\\adaptor.mars.external-1715699429.4551916-471-14-75aa5f3d-f211-47b2-8e50-5a18552bb944.nc"}, 
            {'surface':0,'water':0,'wind_dir':0,'wind_speed':0,'water_x':0,'water_y':0}, 
            datetime.datetime.strptime("2018-01-01-12", "%Y-%m-%d-%H"), 
            "out_raster_path",
            in_coast_polyline=r"C:\GIS\CAA_2023\Data\Script_Datasets\Coastline\EU_coast_1km_select.shp", 
            coast_dist_limit=None)
                