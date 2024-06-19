'''
Created on May 26, 2024

@author: Karl
'''
import cdsapi

c = cdsapi.Client()

c.retrieve(
    'reanalysis-uerra-europe-single-levels',
    {
        'origin': 'mescan_surfex',
        'variable': '10m_wind_speed',
        'year': '2018',
        'month': '01',
        'day': '01',
        'time': [
            '00:00', '06:00', '12:00',
            '18:00',
        ],
        'format': 'netcdf',
    },
    'download.nc')