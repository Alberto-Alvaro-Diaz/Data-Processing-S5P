#TROPOMI S5P OUTPUTS, VARIABLES AND PLOTS

#IMPORTS
#General Modules
import numpy as np

#For plotting S5P soundings 
from matplotlib.patches import Polygon

#scripts
import TropomiPlot as TP

PLOT=True

def tropomi_outputs(File, Lats, Long, region_area , plot=PLOT):
    #Inputs
        #Files: File path of Tropomi S5P data
        #Lats: Latitude coordinates of the points of interest (City, Landfills...)
        #Long: Longitude coordinates of the points of interest 
        #plot: Verification to plot the overpass and the precision 
    
    #corners area of the soundings
    lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
    lon_min, lon_max = region_area['lon_min'], region_area['lon_max']
        
    #Load S5P observations from file
    qa = File.groups['PRODUCT']["qa_value"][:] #quality value 
    time = File.groups['PRODUCT']["time_utc"][:]
    var = File.groups['PRODUCT']["methane_mixing_ratio_bias_corrected"][:]
    varerr = File.groups['PRODUCT']["methane_mixing_ratio_precision"][:]
    lat = File.groups['PRODUCT']["latitude"][:]
    lon = File.groups['PRODUCT']["longitude"][:]
    lat0 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds'][0,:,:,0]
    lat1 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds'][0,:,:,1]
    lat2 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds'][0,:,:,2]
    lat3 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['latitude_bounds'][0,:,:,3]
    lon0 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds'][0,:,:,0]
    lon1 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds'][0,:,:,1]
    lon2 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds'][0,:,:,2]
    lon3 = File.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']['longitude_bounds'][0,:,:,3]
    
    #Time is stored for each scanline ("swath"). Assign time to every pixel.
    #print(np.max(qa))
    time, dummy = np.meshgrid(time[0],qa[0,0,:])
    time = time.T #Due to convention of meshgrid.
    #Other potentially interesting variables, which are available in the file: surface pressure, water column, column averaging kernel
        
    #Create list of soundings
    qa = np.ndarray.flatten(qa)
    time = np.ndarray.flatten(time)
    var = np.ndarray.flatten(var)
    varerr = np.ndarray.flatten(varerr)
    lat, lon = np.ndarray.flatten(lat) , np.ndarray.flatten(lon)
    lat0 , lat1, lat2, lat3 = np.ndarray.flatten(lat0), np.ndarray.flatten(lat1), np.ndarray.flatten(lat2), np.ndarray.flatten(lat3)
    lon0 , lon1, lon2, lon3 = np.ndarray.flatten(lon0), np.ndarray.flatten(lon1), np.ndarray.flatten(lon2), np.ndarray.flatten(lon3)
    
    #Filter the chosen region and time
    quality_product_min=0.1
    mask = [(lat>lat_min)& (lat<lat_max)& (lon_min<lon)&(lon<lon_max)&(qa>quality_product_min)][0]
    time = time[mask]
    var = var[mask]
    varerr = varerr[mask]
    lat, lon = lat[mask] ,lon[mask]
    #cornes of the sounding area
    lat_bounds = np.array([lat0[mask], lat1[mask], lat2[mask], lat3[mask]])
    lon_bounds = np.array([lon0[mask], lon1[mask], lon2[mask], lon3[mask]])

   
    #Find the time of the overflight
    if len(time)>0 or len(var)>0:
        time_Ov = time[int(len(time)/2)]
        print('Time of overflight:',time_Ov,'\n')
        date =f'{time_Ov[0:4]}-{time_Ov[5:7]}-{time_Ov[8:10]}-{time_Ov[11:19]}'        
        
        #Create a patch with the 4 corner points for every sounding.
        patches = []
        for k in range(len(var)):
            poly = np.array([\
                [lon_bounds[0][k],lat_bounds[0][k]],
                [lon_bounds[1][k],lat_bounds[1][k]],
                [lon_bounds[2][k],lat_bounds[2][k]],
                [lon_bounds[3][k],lat_bounds[3][k]]\
                ])
            polygon = Polygon(poly)
            patches.append(polygon)   
    
        #variables to plot
        if plot==True:
            data=[var,varerr] 
            for i,variable in enumerate(data):   
                if i==0:
                    title=f'S5P - CH4 measurements ({date})'
                    colorbar_label='XCH4 [ppb]'
                else:
                    title=f'S5P - CH4 precision measurements ({date})'
                    colorbar_label='XCH4 precision [ppb]'
        
                graph=TP.tropomi_plot(variable, Lats, Long, patches, region_area, title, colorbar_label)  
    else:
        print('There is not data in the region sought')
        patches=1  
        date=0

    return var, lon, lat, lon_bounds, lat_bounds, patches,date
   
