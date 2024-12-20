#DOWNLOA DATA TROPOMI AUTOMATIC QUICK LOOK

#Load modules
import os
import pandas as pd
import requests
import zipfile
import shutil

#General modules
import numpy as np
from netCDF4 import Dataset

#scripts
import TropomiData as TD

import warnings
warnings.filterwarnings("ignore") #ignore warnings

#Here the S5P data is downloaded atomatically for the years, months and days selected.
#For select the overpasses a condition criteria is used with 3 conditions:
    #It exists the point of interest pixel
    #It exists at least a n number of pixel around of our point of interest
    #It exists at least a n number of pixel in the all image

def vectorial_product(p1,p2,p3): #for condition 2, see below
    return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])

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

#Use your credentials for the Copernicus Data Space Ecosystem
username = "username1" 
username2 = "username2"
password = "password"

#Variant with importing credentials from a file:
#from creds import * #The file "creds" contains username and password.

#Input parameters for the request

#Time scale
years=['2022'] #'2022','2023','2024'
months=['02'] #'06','07',
initial_days=['01'] 
final_days=['03'] 
#initial_days=['20']
#final_days=['22']

all_files_saved={}
for year in years:
    for month in months:
        print(year,'-',month,'data will be downloaded:\n')
        for i in range(len(initial_days)):
            try:
                start_date = year+'-'+month+'-'+initial_days[i]
                end_date = year+'-'+month+'-'+final_days[i]
                #Search for products of a specific "data_collection". (E.g. Sentinel-2, Sentinel-5P, ...).
                data_collection = "SENTINEL-5P"
                
                #Set the "contains" argument. (Example: The 'Name' variable should contain 'CH4').
                gas='CH4'
                contains_argument = "contains(Name,'S5P_OFFL_L2__"+ gas +"')"
                
                #Set a polygon for the search (lon1 lat1, lon2 lat2, ..., lon1 lat1)
                aoi = "POLYGON((-5 43,-1 43,-1 37,-5 37,-5 43))'" #pinto 40.257118 Longitud: -3.637550
                
                #Get the catalogue of Sentinel files using an individual request configuration. (Request with the "contains" argument and with coordinates).
                #https://documentation.dataspace.copernicus.eu/APIs/OData.html#filter-option
                json = requests.get(
                    f"https://catalogue.dataspace.copernicus.eu/odata/v1/Products?$filter=Collection/Name eq '{data_collection}' and ContentDate/Start gt {start_date}T00:00:00.000Z and ContentDate/Start lt {end_date}T00:00:00.000Z and {contains_argument} and OData.CSC.Intersects(area=geography'SRID=4326;{aoi})" 
                    ,verify=False
                ).json()
                
                #create the de folders
                Folder_name_final= 'C:/M-MRS/SP5_Data/'+year+'_'+gas+'/'
                if not os.path.exists(Folder_name_final):
                    os.makedirs(Folder_name_final)
                
                #Get the catalogue as a pandas DataFrame
                df = pd.DataFrame.from_dict(json['value'])
                
                #Print the catalogue with a limited number of entries
                
                rows=df.head(30).shape[0]
                columns=df.head(30).shape[1]
                for i in range(rows):
                    print(i, 'Name:' , df['Name'][i],'\n', 
                          'Origin Date:' , df['OriginDate'][i],'\n',
                          'Content Date:' , df['ContentDate'][i],'\n')
                
                #Select an ID from the catalogue
                #id_name = input('Select the ID of the data for download it: ')
                for id_name in range(rows):
                    
                    #Path folder
                    name_data_folder=df['Name'][(int(id_name))][:-3]
                    Folder_name_original= 'C:/M-MRS/scripts/scripts_new/'+name_data_folder
                    
                    #print(df['Footprint'][2])
                    print('File #',id_name,'-------------------------')
                    ID = df['Id'][int(id_name)]
                    
                    if os.path.exists(Folder_name_final+name_data_folder):
                        print(f"\nAlready exists: {name_data_folder}, skipping to the next file\n")
                        
                    else: 
                        #Access to data via APIs requires an access token, which has to be aquired in a first step:
                        
                        #Template to get the access token.
                        def get_access_token(username: str, password: str) -> str:
                            data = {
                                "client_id": "cdse-public",
                                "username": username,
                                "password": password,
                                "grant_type": "password",
                            }
                            try:
                                r = requests.post(
                                    "https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token",
                                    data=data, verify=False
                                )
                                r.raise_for_status()
                            except Exception as e:
                                raise Exception(
                                    f"Access token creation failed. Reponse from the server was: {r.json()}"
                                )
                            return r.json()["access_token"]
                        
                        #Get access token
                        access_token = get_access_token(username, password)
                        
                        #Print access token
                        #print(access_token)
                        
                        #Download the individual file having its ID and an access token (based on a code template from the the guide):
                        
                        url = f"https://zipper.dataspace.copernicus.eu/odata/v1/Products({ID})/$value"
                        headers = {"Authorization": f"Bearer {access_token}"}
                        
                        session = requests.Session()
                        session.headers.update(headers)
                        response = session.get(url, headers=headers, stream=True, verify=False)
                        
                        with open("Download_Product.zip", "wb") as file:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    file.write(chunk)
                        
                        
                        #Note: The download produces a zip file.
                        
                        #Extract the zip file
                        with zipfile.ZipFile('Download_Product.zip', 'r') as zip_ref:
                            zip_ref.extractall('./')
                        
                        #Remove the zip file
                        os.remove("Download_Product.zip")
                        
                        #Note: As a result, the product file is stored within the same directory as this notebook. 
                        #There are actually two files: A directory which contains a cdl file and the product's netcdf file. I don't make use of the cdl file.
                
                        shutil.move(Folder_name_original, Folder_name_final)
                        print(f"\nSaved: {name_data_folder}")
                    
                        #quick look the data
                        
                        File = Dataset(f'{Folder_name_final}/{name_data_folder}/{name_data_folder}.nc')
                        
                        #function
                        [var,lon,lat,long_bounds,lat_bounds, patches, date]=TD.tropomi_outputs(File, Lats, Long, region_area, plot=False)
                        
                        
                        #if there are data
                        if len(var)>0:
                            #check there are enough data
                            #Condition 1: there are at least # pixels surrounding the landfill
                            n_min=75
                            radius=0.5 #deg
                            #center of coordinates our point of interest
                            lat_new , lon_new = lat-Lats[0] , lon-Long[0]
                            mod=np.sqrt(lat_new**2+lon_new**2)
                            cont=0
                            for i in range(len(lat)):
                                if mod[i]<radius:
                                    cont+=1
                                    
                            Cond1=False
                            if cont>n_min:
                                Cond1=True
                                print('Condition 1 satisfied, #pixels surrounding landfill:', cont)
                            else:
                                print('Condition 1 NOT satisfied, #pixels surrounding landfill:', cont, 'min:', n_min)
                            
                            #Condition 2: the pixel landfill exists
                            Cond2=False
                            idx=[]
                            for i in range(len(long_bounds[0])): #the vectorial of the pixels corners and the point of interest is performed in order to know if the point is inside of any pixel
                                v1= vectorial_product([long_bounds[0][i],lat_bounds[0][i]],[long_bounds[1][i],lat_bounds[1][i]],[Long[0],Lats[0]])
                                v2= vectorial_product([long_bounds[1][i],lat_bounds[1][i]],[long_bounds[3][i],lat_bounds[3][i]],[Long[0],Lats[0]])
                                v3= vectorial_product([long_bounds[3][i],lat_bounds[3][i]],[long_bounds[0][i],lat_bounds[0][i]],[Long[0],Lats[0]])
                                if v1>=0 and v2>=0 and v3>=0:
                                    idx.append(i)
                                if v1<=0 and v2<=0 and v3<=0:
                                    idx.append(i)
                                v4= vectorial_product([long_bounds[1][i],lat_bounds[1][i]],[long_bounds[2][i],lat_bounds[2][i]],[Long[0],Lats[0]])
                                v5= vectorial_product([long_bounds[2][i],lat_bounds[2][i]],[long_bounds[3][i],lat_bounds[3][i]],[Long[0],Lats[0]])
                                v6= vectorial_product([long_bounds[3][i],lat_bounds[3][i]],[long_bounds[1][i],lat_bounds[1][i]],[Long[0],Lats[0]])
                                if v4>=0 and v5>=0 and v6>=0:
                                    idx.append(i)
                                if v4<=0 and v5<=0 and v6<=0:
                                    idx.append(i)
                            if idx: 
                                Cond2=True 
                                print('Condition 2 satisfied, landfill pixel exists')
                            else:
                                print('Condition 2 NOT satisfied, landfill pixel does not exist')
                                
                            #Condition 3: there area at leat # pixels in the total image
                            Cond3=False
                            n_min=170
                            if len(var)>n_min:
                                Cond3=True
                                print('Condition 3 satisfied, #pixels in the image:', len(var))
                            else:
                                print('Condition 3 NOT satisfied, #pixels in the image:', len(var), 'min:', n_min)
                            
                            #all conditions must be satisfied
                            if Cond1==False or Cond2==False or Cond3==False:
                                print('Conditions did not pass the filters')
                                File.close()
                                os.remove(Folder_name_final+name_data_folder+'/'+name_data_folder+'.nc')
                                shutil.rmtree(Folder_name_final+name_data_folder)
                                print(f"Removed folder\n")
                            else:
                                print('Conditions did pass the filters, file saved!\n')
                                all_files_saved[date] = name_data_folder
                                
                            
                        else:
                            print(f"\n{Folder_name_final+name_data_folder} There is not data")
                            File.close()
                            os.remove(Folder_name_final+name_data_folder+'/'+name_data_folder+'.nc')
                            shutil.rmtree(Folder_name_final+name_data_folder)
                            print(f"Removed folder\n")
            except Exception as e:
                print(f"Error raised {e}, follow with the next username")
                username=username2
                    
            
print('Files saved: ',all_files_saved,'\n')
print('Final usearname used:',username,'\n')
print('\nEND')



        
        
    

