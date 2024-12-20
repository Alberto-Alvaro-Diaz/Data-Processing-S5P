#WIND FIELDS AT DIFFERENT PRESSURE LEVELS, ERA5 (PROFILE)
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels?tab=form

#IMPORTS
#General Modules
import numpy as np
import matplotlib.pyplot as plt

#data format
import DateData 

#Plot geo maps
import cartopy.crs as ccrs
import cartopy.feature as cf
from shapely.geometry import box
from cartopy.feature import ShapelyFeature
from matplotlib.patches import Patch

PATH_SAVE=[]

def vertical_profiel_wind_field(File,Lats,Long,path_save1=PATH_SAVE,path_save2=PATH_SAVE):
    #Inputs
        #File: File path of data from ERA5
        #Lats: Latitude coordinates of the points of interest (City, Landfills...)
        #Long: Longitude coordinates of the points of interest
        #path_save: path and figure name where save the plot
              
    #calculation of data date 
    date=DateData.date_data(File)
    print('Date of ERA5 Profile data ',date,'\n')
    
    #latitudes and longitudes array
    lat_array, long_array=np.array(File['latitude']) ,np.array(File['longitude'])
    
    #calculate wind velocity in area surrounding the landfill
    array_x=[]
    array_y=[]
    lats=lat_array-Lats[0]
    long=long_array-Long[0]
    radius= 0.5 #degrees coordinates
    for i in range(len(lat_array)):
        for j in range(len(long_array)):
            if np.sqrt(lats[i]**2+long[j]**2) <= radius:
                array_x.append(j)
                array_y.append(i)
    
    #coordinates Pinto's landfill matrix units
    point_y, point_x= np.argmin(abs(lat_array-Lats[0])),  np.argmin(abs(long_array-Long[0]))
    
    #wind velocities profile
    #mean of the surrounding points
    # wind_u , wind_v = np.mean(np.array(File['u'][:,:,array_y,array_x]),axis=3), np.mean(np.array(File['v'][:,:,array_y,array_x]),axis=3)
    # wind_u , wind_v = np.mean(wind_u,axis=2), np.mean(wind_v,axis=2)
    
    #mean of the nearest point (valid_time, pressure_level, latitude, longitude)
    wind_u , wind_v = np.array(File['u'][:,:,point_y,point_x]), np.array(File['v'][:,:,point_y,point_x])
    wind_u , wind_v = np.mean(wind_u, axis=0), np.mean(wind_v, axis=0) #temporal mean
    
    #calculate the wind direction for each layer
    angle = np.degrees(np.arctan(wind_v/wind_u)) #angle 0 set to positive u direction (+ x direction)
     
    #calculate wind module velocities at each layer
    wind_module=np.sqrt(wind_u**2+wind_v**2)
    
    #layers level and calculate height
    layers_pressure=np.array(File['pressure_level'][:])
    # rho=1.204 #kg/m^3 assumed cte
    # h=100*layers_pressure/(rho*g)
    # h=abs(h-h[0])  #set sea level
    
    #temperature profile
    temperature=np.array(File['t'][:,:,point_y,point_x])
    temperature=np.mean(temperature, axis=0) #temporal mean
    
    #humidity
    humidity=np.array(File['r'][:,:,point_y,point_x])
    humidity=np.mean(humidity, axis=0)
    
    #ozone
    ozone=np.array(File['o3'][:,:,point_y,point_x])
    ozone=np.mean(ozone, axis=0)
    
    #hight profile
    M=0.0289644 #kg/mol
    R=8.3144598 #J/(mol*K)
    g=9.81 #gravity m/s^2
    T=temperature[0]
    
    # print(T,layers_pressure[0])
    # print(layers_pressure)
    h=-R/(g*M)*T*np.log(layers_pressure/layers_pressure[0])
    # print(h)
    
    #plot temperature, ozone and humidity profilee
    fig, ax=plt.subplots(figsize=[15,12])  
    
    #ax
    ax.plot(ozone, h,'o-',color='b',markersize=10,label='Ozone Profile')
    ax.set_xlabel('Mass Mixing Ratio (XO3)',fontsize=20, color='b')
    ax.set_ylabel('Altitude (m)',fontsize=20)
    
    ax1=ax.twiny()
    ax1.plot(humidity,h,'o-',color='r',markersize=10,label='Humidity Profile')
    ax1.set_xlabel('Relative Humidity (%)',fontsize=20, labelpad=10, color='r')
    
    ax2=ax.twiny()
    ax2.plot(temperature,h,'o-',color='k',markersize=10,label='Temperature Profile')
    ax2.set_xlabel('Temperature (K)',fontsize=20, labelpad=35, color='k')
    
    
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(axis='x', labelsize=18, labelcolor='b')
    
    ax1.tick_params(axis='x', labelsize=18,labelcolor='r')
    ax2.tick_params(axis='x', labelsize=18,labelcolor='k')

    lines_0, labels_0 = ax.get_legend_handles_labels()
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_0 + lines_1 + lines_2, labels_0 + labels_1 + labels_2, loc='upper center', fontsize=15)
    
    
    plt.title(f'ERA5 - Temperature, Humidity and Ozone Profile ({date})',fontsize=20)
    #labels, colors and sizes
    ax.tick_params(axis='y', labelsize=18)
    ax.tick_params(axis='x', labelsize=18, labelcolor='b')
    
    #grids
    ax.grid(color='grey')
    if path_save2:
        plt.savefig(path_save2)
    plt.show()

    for i in range(2):
        #plot
        fig, ax=plt.subplots(figsize=[15,12])  
        
        #ax
        ax.plot(wind_module, layers_pressure,'*-',color='y')
        plt.gca().invert_yaxis()
        ax.set_xlabel('Wind Velocity (m/s)',fontsize=20, color='b')
        ax.set_ylabel('Pressure level (hPa)',fontsize=20)
        
        #ax1
        ax1=ax.twinx()
        ax1.plot(wind_module, h,'o-',markersize=10, label='Wind Profile', color='b')
        ax1.set_ylabel('Height (m)',fontsize=20)
        
        #ax2
        ax2=ax.twiny()
        ax2.plot(angle,layers_pressure,'o-',markersize=10, label='Wind Direction Profile', color='g')
        ax2.set_xlabel("Wind Direction (º)", fontsize=20, color='g')
        
        plt.title(f'ERA5 - Winds Profile ({date})',fontsize=20)
        #labels, colors and sizes
        ax.tick_params(axis='y', labelsize=18)
        ax.tick_params(axis='x', labelsize=18, labelcolor='b')
        ax1.tick_params(axis='y', labelsize=18)
        ax2.tick_params(axis='x', labelsize=18,labelcolor='g')
        
        #grids
        ax.grid(color='k',axis='x')
        ax1.grid(color='k')
        ax2.grid(color='g')
          
        if i==1:
            #integrate til some height
            N=n #layer number (points in the graph) up to 1500 meters
            ax1.plot(wind_module[:N],h[:N],'o-',markersize=10, color='r', label='Integration points')
        
        #legends from both x axis en ona
        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper center', fontsize=15)
        
        if path_save1:
            plt.savefig(path_save1)

        plt.show()
        
        if i==0:
            n=int(input('Select PBL altitude (# layer): '))
            #n=3
            PBL_altitude=round(h[n-1],3)
            print('Altitude selected:',PBL_altitude,'m\n')
    

    #calculate the mean wind velocity and errors
    std_wind_u , std_wind_v = np.std(wind_u[:N]) , np.std(wind_v[:N]) #std temporal profile 
    
    wind_u, wind_v =np.mean(wind_u[:N]), np.mean(wind_v[:N])
    print(f'Wind Field Mean Vectors up to {PBL_altitude}m x,y: ',wind_u,wind_v)    
    wind_mean=np.sqrt(wind_u**2+wind_v**2)
    print(f'Wind Mean Module up to {PBL_altitude}m', wind_mean,'\n')
    
    delta_U = abs(wind_u/wind_mean) *std_wind_u + abs(wind_v/wind_mean) * std_wind_v
    
    return wind_u , wind_v ,  wind_mean , PBL_altitude , delta_U #np.mean(wind_module[:N])


if __name__=='__main__':
    
    from netCDF4 import Dataset
    Lats=np.array([40.257118, 40.336427, 40.457619, 39.854100, 40.41831])
    Long=np.array([-3.637550, -3.590375, -3.363080, -4.168100, -3.70275])
    
    Folder_name_profile='C:/M-MRS/ERA5_Data/Wind_profile/2024/'
    File_name='Pressure_Levels_2024_06_01'
    File_name_profile=Folder_name_profile+File_name
    
    File= Dataset(f'{File_name_profile}.nc')
    
    (wind_u , wind_v ,  wind_mean , PBL_altitude , delta_U)=vertical_profiel_wind_field(File,Lats,Long)
    print(wind_u , wind_v ,  wind_mean , PBL_altitude , delta_U)
    
    
    
    
    
    
    
    
    
    