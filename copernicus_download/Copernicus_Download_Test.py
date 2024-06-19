'''
Created on May 10, 2024

@author: Karl
'''
import arcpy
import copernicusmarine as cm

def copernicus_download(extent_fc, start_datestring, end_datestring, out_name):
    
    extent_desc = arcpy.Describe(extent_fc)
    extent_obj = extent_desc.extent
    
    # North Sea sea surface height.
    #cmems_mod_nws_phy-ssh_my_7km-2D_PT1H-i
    
    # North Sea currents
    
    # World Dataset
    #cmems_mod_glo_phy_anfc_0.083deg_PT1H-m
    
    #u-v dataset details
    #dataset_id="cmems_mod_nws_phy-uv_my_7km-2D_PT1H-i",
    #  variables=["uo", "vo"],
    
    # Download dataset.
    cm.subset(
      dataset_id="cmems_mod_nws_phy-ssh_my_7km-2D_PT1H-i",
      variables=["zos"],
      minimum_longitude=extent_obj.XMin,
      maximum_longitude=extent_obj.XMax,
      minimum_latitude=extent_obj.YMin,
      maximum_latitude=extent_obj.YMax,
      start_datetime=start_datestring,
      end_datetime=end_datestring,
      minimum_depth=0,
      maximum_depth=1,
      output_filename = out_name,
      output_directory = r"C:\GIS\Maritime_Encounters\OrmeSim\NetCDF"
    )
    
# Sample datestring: "2018-01-01T00:00:00"
# Get extent from shapefile (must be projected in WGS).
extent_fc = r"C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Features\DEM_footprint_WGS84"
copernicus_download(extent_fc, "2018-01-01T00:00:00", "2018-12-31T23:59:59", "current_zos_2018.nc")