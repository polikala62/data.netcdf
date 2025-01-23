'''
Created on 20 Nov 2024

@author: Karl Smith
'''

import arcpy
arcpy.env.overwriteOutput = True

points = r"C:\GIS\ArcPro_Projects\PM_Med_Survey\Temp\test_obs_pts_01.shp"
dem = r"C:\GIS\ArcPro_Projects\PM_Med_Survey\Rasters\gmted_mea150_clip_02.tif"
vshed = r"C:\GIS\ArcPro_Projects\PM_Med_Survey\Temp\tempvs.tif"

# Create update cursor for feature class 
with arcpy.da.UpdateCursor(points, ["SHAPE@", "coast_vis"]) as cursor:
    for row in cursor:
        # Do viewshed.
        arcpy.ddd.Viewshed2(
        in_raster=dem,
        in_observer_features=row[0],
        out_raster=vshed,
        analysis_type="FREQUENCY",
        refractivity_coefficient=0.13,
        surface_offset="0 Meters",
        observer_elevation="5 Meters",
        analysis_method="PERIMETER_SIGHTLINES",
        analysis_target_device="GPU_THEN_CPU"
        )
        
        # Convert viewshed to a numpy array
        array = arcpy.RasterToNumPyArray(vshed, nodata_to_value=0)

        # Sum the array
        vshed_sum = int(array.sum())
        print(vshed_sum)
        # Delete raster.
        #arcpy.management.Delete(vshed)
        
        # Update row.
        row[1] = vshed_sum
        cursor.updateRow(row)
        
print("Script finished!")