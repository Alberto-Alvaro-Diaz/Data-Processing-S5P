#WIND FIELDS AT 10 AND 100m FROM THE GROUND AND SURFACE PRESSURE 
#https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form

#IMPORTS
#General Modules
import numpy as np
import matplotlib.pyplot as plt
import re

#Plot geo maps
import cartopy.crs as ccrs
import cartopy.feature as cf

#scripts
import DateData  

def surface_pressure_wind_field(File,Lats,Long,region_area): 
    #Inputs
        #Files: File path of data from ERA5 (SP,10,100)
        #Lats: Latitude coordinates of the points of interest (City, Landfills...)
        #Long: Longitude coordinates of the points of interest 
        #region_area: Dictionary which defines the limit of the area (cornes latidue and longitudes)
    
    #corners area of the soundings
    lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
    lon_min, lon_max = region_area['lon_min'], region_area['lon_max']
        
    variables_dictionary={'u10':np.mean(File['u10'][:], axis=0),'u100':np.mean(File['u100'][:], axis=0),
                          'v10':np.mean(File['v10'][:], axis=0),'v100':np.mean(File['v100'][:], axis=0),
                          'sp':np.mean(File['sp'][:], axis=0)}
       
    #variables definition
    lat, lon =np.array(File['latitude'][:]) , np.array(File['longitude'][:])
    
    #the array of Lats and Long is compared with the Landfills facilities coordinates (Pinto in this case)
    point_y, point_x = np.argmin(abs(lat-Lats[0])),np.argmin(abs(lon-Long[0]))  #The coordinates of Pinto's facility is obtained in matrix indices units
    
    #corners plot in matrix indices
    lat_min_indx, lat_max_indx, lon_min_indx, lon_max_indx = np.argmin(abs(lat-lat_min)), np.argmin(abs(lat-lat_max)), np.argmin(abs(lon-lon_min)), np.argmin(abs(lon-lon_max))

    date=DateData.date_data(File)
    print('Date of ERA5 Surface Pressure, Wind 10m and 100m ',date)
                
    #surface pressure
    surface_pressure=File['sp'][0]
    print('Surface Presure (Pa) (Pinto): ', surface_pressure[point_x,point_y],'\n')
    
    #PLOTS
    height=[] #'10','100'
    wind_fields={}
    for i in range(len(height)+1):   
        fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[12,14])
        if i==0:
            #variables definition
            surface_pressure_plot=surface_pressure[lat_max_indx:lat_min_indx+1,lon_min_indx:lon_max_indx+1]
            [X, Y] = np.meshgrid(lon[lon_min_indx:lon_max_indx+1], lat[lat_max_indx:lat_min_indx+1])
            
            #data      
            plt1 = ax.contourf(X, Y, surface_pressure_plot/100, cmap='nipy_spectral_r', transform=ccrs.PlateCarree())
           
            #colorbar
            cbar_ax = fig.add_axes([0, 0, 0.1, 0.1]) # Dummy values prior to finetuning the cbar position
            pos = ax.get_position() # Get the axes position
            cbar_ax.set_position([pos.x0 + pos.width + 0.01, pos.y0, 0.04, pos.height])
            cbar = plt.colorbar(plt1, cax=cbar_ax)
            cbar_ax.tick_params(labelsize=8)
                  
            ax.set_title(f'ERA5 - Surface Pressure (hPa) ({date})')         
            
        else:
            #variables definition
            u=np.array(variables_dictionary['u'+height[i-1]]) #eastward component of the wind (m/s)
            v=np.array(variables_dictionary['v'+height[i-1]]) #northward component of the wind (m/s)
                       
            #data
            ax.barbs(lon[:],lat[:], u[:,:], v[:,:],transform=ccrs.PlateCarree())
            
            plt.title(f'ERA5 - Winds Fields {height[i-1]} meters ({date})')
            
        
        #plot features    
        ax.add_feature(cf.COASTLINE,linewidth=1,zorder=4)
        ax.add_feature(cf.BORDERS,linewidth=1,zorder=3)
        ax.add_feature(cf.STATES, linewidth=1, edgecolor='grey',zorder=3)
        ax.set_extent([lon_min, lon_max, lat_min, lat_max])
        
        #grid features 
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=1,linestyle='--', color='gray', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {'size': 12}
        gl.ylabel_style = {'size': 12}
        
        #Madrid City and Facilities
        ax.plot(Long[:-1],Lats[:-1],'k*', markersize=5, transform=ccrs.Geodetic(), label='Landfill Facilities')
        ax.plot(Long[-1],Lats[-1],'r*', markersize=5, transform=ccrs.Geodetic(),label='Madrid City')
        
        ax.legend(loc='upper left')
        plt.show()
    
    return surface_pressure[point_x,point_y], [variables_dictionary['u10'],variables_dictionary['u100'], variables_dictionary['v10'], variables_dictionary['v100']], lat, lon







