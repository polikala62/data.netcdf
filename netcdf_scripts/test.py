'''
Created on Mar 21, 2024

@author: Karl
'''
import netCDF4 as nc
import numpy as np
import xarray as xr

def check_netcdf(ncdf_path):
    
    # Create dataset.
    rootgrp = nc.Dataset(ncdf_path, "r", format="NETCDF4")
    variables = rootgrp.variables
    dimensions = rootgrp.dimensions
    
    print(dimensions)
    print()
    print(variables)
    
    '''
    sub_area = variables['sub_area']
    print(sub_area.dimensions)
    print(sub_area.shape)
    
    print(np.average(sub_area,axis=0).shape)
    #print(rootgrp.dimensions)
    #print(rootgrp['z'][:])
    #variables = rootgrp.variables['sub_area'][:][0]
    print()
    xvar = variables['easting']
    print(np.sum(xvar))
    print(sub_area.size)
    
    #print(variables)
    '''
check_netcdf(r"C:\GIS\ArcPro_Projects\Visibility_Testing\Output\med_test_01.nc")
#check_netcdf(r"C:\GIS\Maritime_Encounters\OrmeSim\Output\full_area_03.nc")
#check_netcdf(r"C:\GIS\Maritime_Encounters\OrmeSim\NetCDF\Orme_Water_2018_01_01.nc")
#check_netcdf("C:\\GIS\\Maritime_Encounters\\OrmeSim\\Output\\full_area_03.nc")
#check_netcdf(r"C:\GIS\Maritime_Encounters\Copernicus_Testing\Download_Test\Download\adaptor.mars.external-1715699429.4551916-471-14-75aa5f3d-f211-47b2-8e50-5a18552bb944.nc")