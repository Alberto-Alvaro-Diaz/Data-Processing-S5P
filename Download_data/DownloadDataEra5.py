#DOWNLOAD DATA FORM ERA5

#imports
#General modules
import numpy as np
from netCDF4 import Dataset

#API files and folders
import cdsapi
import shutil
import os

#scripts
import TropomiData as TD

import warnings
warnings.filterwarnings("ignore") #ignore warnings

#Here the data of ERA5 is downloaded automatically. This script uses the folder wich were saved the
#files of S5P to extract the date of the overpasess to download the data.
#Here are downloaded two files:
    #Single level (pressure, 10 and 100 meters wind fields)
    #Pressure levels for wind fields and temperature.

####################

#Select region (Pinto-Madrid-Toledo) corners area of the soundings
region_area= {'lat_min':39.75,'lat_max':41.25, 'lon_min':-4.75,'lon_max':-3 }
lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
lon_min, lon_max = region_area['lon_min'], region_area['lon_max']

#Facilities coordinates and Madrid city: According to PRTR España (web)
#Pinto coordinates: Latitud: 40.257118 Longitud: -3.637550
#Valdemingomez coordinates: Latitud 40.336427 Longitud: -3.590375
#Acalá coordinates Latitud: Latitud 40.457619 Longitud: -3.363080 
#Toledo ESTACIÓN DE TRATAMIENTO DE RESIDUOS Latitud: 39,854100 Longitud: -4,168100
#Madrid city coordinates: 40.41831, -3.70275
Lats=np.array([40.257118, 40.336427, 40.457619, 39.854100, 40.41831])
Long=np.array([-3.637550, -3.590375, -3.363080, -4.168100, -3.70275])

####################
years=['2024'] #['2024','2023','2022','2021']
for year in years:
    folder_name=f"C:/M-MRS/SP5_Data/{year}_CH4" #where the tropomi data is downloaded
    folder_name_downloaded="C:/M-MRS/scripts/scripts_new/" #where the ERA5 will be downloaded 
    folder_name_single="C:/M-MRS/ERA5_Data/Surface_pressure_and_10_100_m_wind_comparision" #where we want to save ERA5 data single level
    folder_name_pressure="C:/M-MRS/ERA5_Data/Wind_profile_comparision" #where we want to save ERA5 data Pressure levels
    files=os.listdir(folder_name)
    dates=[]
    hours=[]
    
    for file in files:
        File = Dataset(f'{folder_name}/{file}/{file}.nc')
        [var,lon,lat,long_bounds,lat_bounds, patches, date]=TD.tropomi_outputs(File, Lats, Long, region_area, plot=False)
        dates.append(date[:10])
        hours.append(date[11:])
    
    #SINGLE LEVEL
    c = cdsapi.Client(verify=False)
    for idx, date in enumerate(dates):
            
        year=date[:4]
        month=date[5:7]
        day=date[8:10]
        time=hours[idx][:2]
        minutes=hours[idx][3:5]
        if int(minutes) >=30:
            time=int(time)+1
        time=str(time)+':00'
        
        time1=str(int(time[:2])-1)+':00'
        time2=str(int(time[:2])-2)+':00'
        time3=str(int(time[:2])-3)+':00'
        
        if not os.path.exists(folder_name_single+'/'+year):
            os.makedirs(folder_name_single+'/'+year)
        file_name_single='Single_Level_'+year+'_'+month+'_'+day+'.nc'  
        file = folder_name_downloaded+'/'+file_name_single
        folder_saved =folder_name_single+'/'+year
        if os.path.exists(folder_saved+'/'+file_name_single):
            print(f"\nAlready exists: {file}")
        else:
            try:
                
                c.retrieve(
                    'reanalysis-era5-single-levels',
                    {
                        'product_type': 'reanalysis',
                        'format':'netcdf',
                        'year': year,
                        'month': month,
                        'day': day,
                        'time': [time,time1,time2,time3],
                        'variable': [
                            'surface_pressure','boundary_layer_height',
                        ],
                        'area': [
                            42, -5, 39, 3,
                        ],
                    },
                    file_name_single
                  )
                
                #Move file to folder
                shutil.move(file, folder_saved)
                print('\nDownloaded file:', file_name_single,'in:',folder_saved,'\n')
                    
            except Exception as e:
                print(f"\n'{e}' with date {year}-{month}-{day}. Follow with the next data")
            
      
    #PRESSURE LEVELS   
    c = cdsapi.Client(verify=False) 
    for idx, date in enumerate(dates):
        year=date[:4]
        month=date[5:7]
        day=date[8:10]
        time=hours[idx][:2]
        minutes=hours[idx][3:5]
        if int(minutes) >=30:
            time=int(time)+1
        time=str(time)+':00'
        
        time1=str(int(time[:2])-1)+':00'
        time2=str(int(time[:2])-2)+':00'
        time3=str(int(time[:2])-3)+':00'
        
        if not os.path.exists(folder_name_pressure+'/'+year):
            os.makedirs(folder_name_pressure+'/'+year)
            
            
        file_name_pressure='Pressure_Levels_'+year+'_'+month+'_'+day+'.nc'    
        file = folder_name_downloaded+'/'+file_name_pressure
        folder_saved =folder_name_pressure+'/'+year  
        if os.path.exists(folder_saved+'/'+file_name_pressure):
            print(f"\nAlready exists: {file}")
        else:
            try:
                #relative_humidity
                c.retrieve(
                    'reanalysis-era5-pressure-levels',
                    {
                        'product_type': 'reanalysis',
                        'format': 'netcdf',
                        'variable': [
                            'u_component_of_wind', 'v_component_of_wind','temperature'
                            ,'relative_humidity','ozone_mass_mixing_ratio', 
                        ],
                        'pressure_level': [
                            '5', '20', '50', '175',
                            '250','300','350',
                            '400','450','500',
                            '550','600','650',
                            '700', '750', '775',
                            '800', '825', '850',
                            '875', '900', '925',
                            '950', '975', '1000',
                        ],
                        'year': year,
                        'month': month,
                        'day': day,
                        'time': [time,time1,time2,time3],
                        'area': [
                            42, -5, 39, 3,
                        ],
                    },
                    file_name_pressure
                    )            
                    
                #Move file to folder
                shutil.move(file, folder_saved)
                print('\nDownloaded file:', file_name_pressure,'in:',folder_saved,'\n')
                
            except Exception as e:
                print(f"\n'{e}' with date {year}-{month}-{day}. Follow with the next data")
            
print('\nEND')
        

