'''
Created on 23 Oct 2024

@author: Karl Smith
'''

import arcpy, math
import numpy as np
from netcdf_scripts import numpy_array, netcdf_to_raster

wgs = arcpy.SpatialReference("WGS 1984")

wind_netcdf = r"C:\GIS\NetCDF_Manipulation\Data\20240101_2dh-CMCC--RFVL-MFSeas8-MEDATL-b20240116_an-sv09.00.nc"
wind_slice_dict = {'time':[0,24], 'lat':[0,380], 'lon':[0,1287]}
#wind_slice_dict = {'time':[0,24], 'lat':[50,250], 'lon':[500,550]}

wind_u_slice = numpy_array.slice_netcdf(wind_netcdf, 'uo', wind_slice_dict)
wind_v_slice = numpy_array.slice_netcdf(wind_netcdf, 'vo', wind_slice_dict)

wind_flat_u_slice = numpy_array.flatten_slice(wind_netcdf, 'uo', wind_u_slice, "STANDARD_DEVIATION", ['lon', 'lat'])
wind_flat_v_slice = numpy_array.flatten_slice(wind_netcdf, 'vo', wind_v_slice, "STANDARD_DEVIATION", ['lon', 'lat'])

wind_velocity_slice = np.sqrt((wind_flat_u_slice**2 + wind_flat_v_slice**2))

netcdf_to_raster.netcdf_to_singleband_raster2(wind_netcdf, 'lon', 'lat', wind_slice_dict, wind_velocity_slice, r'C:\GIS\NetCDF_Manipulation\Output\test_wind_v_01.tif', wgs)

'''
# Read wind u and v variables to arrays.
slice_dict = {'time':[0,24], 'mlev':[0,380], 'lat':[0,96], 'lon':[0,1287]}
pr_slice = numpy_array.slice_netcdf(in_netcdf, 'xl', slice_dict)
slice_avg = numpy_array.flatten_slice_by_dict(in_netcdf, 'xl', pr_slice, {'time':"MEAN", 'mlev':"MEAN"})
'''
