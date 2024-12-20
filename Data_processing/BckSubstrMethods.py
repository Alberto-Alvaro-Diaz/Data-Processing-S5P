#Background substraction different methods
#IMPORTS
#General Modules
import numpy as np

#the following imports are useful for debugs
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cf
from matplotlib.collections import PatchCollection

#scripts
import TropomiFindPlumePixels as TFPP


class  BackSubstractionMethods():
    #Class where is defined the methods to substract the background of the overpass
       
    def rotation_axis(lat, lon, Lats, Long):
        
        #take landfill nearest pixel and find the pixel above or below it. That pixel must be the nearest to the pixel 
        x, y =Long[0]-lon, Lats[0]-lat
    
        idx_x, idx_y=np.argmin(abs(x)), np.argmin(abs(y))
    
        delta_y, delta_x=lat-lat[idx_y], lon-lon[idx_x]
    
        mod=np.sqrt(delta_x**2+delta_y**2)
        value1=np.sort(mod)[0]
        index_landfill = np.where(mod == value1)[0]
        
        mod=np.sqrt((lat-lat[index_landfill])**2+(lon-lon[index_landfill])**2)
        value2=np.sort(mod)[1]
        value3=np.sort(mod)[2]
        index_nearest_landfill=np.where(mod==value2)[0]
        index_nearest_landfill2=np.where(mod==value3)[0]
        
        #find the angle rotation, why two differents angles?
        angle=np.arctan(np.abs(lon[index_nearest_landfill]-lon[index_landfill])/np.abs(lat[index_nearest_landfill]-lat[index_landfill]))  #for long
        angle1=np.pi/2-np.arctan(np.abs(lon[index_landfill]-lon[index_landfill+1])/np.abs(lat[index_landfill]-lat[index_landfill+1]))
        
        #rotate the image to select the column and rows of the pixel above-below and right-left
        lat_new, lon_new=np.zeros_like(lat), np.zeros_like(lon)
        lon_origin, lat_origin=lon-lon[index_landfill], lat-lat[index_landfill]  #select as origin of coordinates the pixel found
        
        #lets apply the rotation
        for i in range(len(lat)):
            lon_new[i]=lon_origin[i]*np.cos(angle)+lat_origin[i]*np.sin(angle)
            lat_new[i]=-lon_origin[i]*np.sin(angle1)+lat_origin[i]*np.cos(angle1)
        #come back to the original coordinates
        lon_new,lat_new =lon_new+lon[index_landfill], lat_new+lat[index_landfill]
        
        return lon_new, lat_new, index_landfill, index_nearest_landfill, index_nearest_landfill2
    
    
    def striping(var, lat, lon, Lats, Long, region_area ,patches, long_bounds, lat_bounds, f):
        #Inputs
            #Lats: Latitude coordinates of the points of interest (City, Landfills...)
            #Long: Longitude coordinates of the points of interest 
            #region area: Dictionary which defines the limit of the area (cornes latidue and longitudes)
            #patches: The pixels of the overpass
            
        #substract background in columns with the aim to avoid striping features
        [lon_new, lat_new, index_landfill, index_nearest_landfill,index_nearest_landfill2]=BackSubstractionMethods.rotation_axis(lat, lon, Lats, Long)
        coords_plume_striping=[]

        #calculate precission for find rows or columns 
        dist_x1, dist_y1 = abs(lon_new[index_landfill+1]-lon_new[index_landfill]) , abs(lat_new[index_landfill]-lat_new[index_nearest_landfill])
        dist_x2, dist_y2 = abs(lon_new[index_landfill-1]-lon_new[index_landfill]) , abs(lat_new[index_landfill]-lat_new[index_nearest_landfill2])
        dist_x ,dist_y = np.mean([dist_x1,dist_x2]), np.mean([dist_y1,dist_y2]) #calculate the minimal distance between the follow and previous landfill pixel
        
        #in order to avoid the issue if one of those pixels do not exist
        if dist_x==0:
            dist_x=0.1 
        if dist_y==0:
            dist_y=0.1
        preccision_x,  preccision_y=np.abs(dist_x)/2 , np.abs(dist_y)/2 #precission x has chosen with the pixel origin and the following 
        
        n_old, n_new ,n = 0, np.pi, 0
        while n_new != n_old and n<30: #limit n to avoid infinite loops
            
            n_old=n_new
            n+=1
            
            bck_pixels=np.empty((2),dtype=object)
            #find the pixels wich have the same long/lat coordinate than the landfill pixel +- precission
            values_x , values_y = np.isclose(lon_new,lon_new[index_landfill],atol=preccision_x/2) , np.isclose(lat_new,lat_new[index_landfill],atol=preccision_y/2)  #columns , rows 
          
            #pixel indices, introduce the arrays in the matrix
            bck_pixels[0]=np.sort(np.array(np.where(values_x)[0])) #column indices
            bck_pixels[1]=np.sort(np.array(np.where(values_y)[0])) #row
            
            all_columns=[]
            for i , i_pixel in enumerate(bck_pixels[0]): #for each column index
                value_y = np.isclose(lat_new,lat_new[i_pixel],atol=preccision_y)
                i_row = np.sort(np.array(np.where(value_y)[0]))  #the row for each column pixel is obained
                
                """
                lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
                lon_min, lon_max = region_area['lon_min'], region_area['lon_max']     
                fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[8,10]  )
                ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
                p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
                p.set_array(var)
                p.set_cmap('coolwarm')
                ax.add_collection(p)
                ax.add_feature(cf.STATES, linewidth=0.2, edgecolor='grey',zorder=3)
                plt.colorbar(p,ax=ax,label='caca',orientation='horizontal',shrink=0.6,pad=0.04)
                ax.plot(lon[i_row],lat[i_row],'r*', markersize=5)                 
                plt.show()
                """
                
                for j, j_pixel in enumerate(i_row): 
                    value_x = np.isclose(lon_new,lon_new[j_pixel],atol=preccision_x/2)
                    j_column = np.sort(np.array(np.where(value_x)[0]))
                    all_columns.append(j_column)
                    
                    """
                    lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
                    lon_min, lon_max = region_area['lon_min'], region_area['lon_max']    
                    fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[8,10]  )
                    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
                    p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
                    p.set_array(var)
                    p.set_cmap('coolwarm')
                    ax.add_collection(p)
                    ax.add_feature(cf.STATES, linewidth=0.2, edgecolor='grey',zorder=3)
                    plt.colorbar(p,ax=ax,label='caca',orientation='horizontal',shrink=0.6,pad=0.04)
                    ax.plot(lon[i_row],lat[i_row],'g*', markersize=5)
                    ax.plot(lon[j_column],lat[j_column],'r*', markersize=5)                  
                    plt.show()
                    """
                    
            columns_sub=list(map(list, dict.fromkeys(tuple(arr) for arr in all_columns)))   
            columns_sub1=np.empty(len(columns_sub),dtype=object)
            for i in range(len(columns_sub1)): 
                columns_sub1[i]=np.array([idx for idx in columns_sub[i] if idx not in coords_plume_striping]) #delete the plume pixels for calculate the background
                
            var_bck_corrected_striping=np.zeros_like(var)
            for idx, column in enumerate(columns_sub):
                bck=np.mean(var[columns_sub1[idx]])  #[-n:]#def background for each column/stripe
                #print(bck) #the bck for each column
                var_bck_corrected_striping[column]=var[column]-bck  #subtract bck to the image
            
            #print(np.mean(var_bck_corrected_striping[coords_plume_striping]))
            coords_plume_striping , bckg_mean = TFPP.find_plume_pixels(var_bck_corrected_striping,lon,lat,Long,Lats,long_bounds, lat_bounds,F=f)

            n_new=len(coords_plume_striping)

        """
        lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
        lon_min, lon_max = region_area['lon_min'], region_area['lon_max']
                
        fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[8,10]  )
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
        p.set_array(var)
        p.set_cmap('coolwarm')
        ax.add_collection(p)
        ax.add_feature(cf.STATES, linewidth=0.2, edgecolor='grey',zorder=3)
        plt.colorbar(p,ax=ax,label='caca',orientation='horizontal',shrink=0.6,pad=0.04)
        
        ax.plot(lon[bck_pixels[0]],lat[bck_pixels[0]],'r*', markersize=5)
        ax.plot(lon[bck_pixels[1]],lat[bck_pixels[1]],'r*', markersize=5)
               
        #colors = plt.cm.rainbow(np.linspace(0, 1, len(columns_sub)))
        #for i, column in enumerate(columns_sub1):
            #ax.plot(lon[column],lat[column], color=colors[i], marker='*', markersize=5)
            
        ax.plot(lon[index_landfill],lat[index_landfill],'b*', markersize=5)
        ax.plot(lon[index_nearest_landfill],lat[index_nearest_landfill],'k*', markersize=5)
        ax.plot(lon[index_nearest_landfill2],lat[index_nearest_landfill2],'k*', markersize=5)
        
        #ax.plot(lon[coords_plume_striping],lat[coords_plume_striping],'bx', markersize=9)
             
        plt.show()
        """
        
        return var_bck_corrected_striping
             
    
        
    def up_wind(var, lat, lon, Lats, Long, region_area ,patches, coords_plume):
        #Inputs
            #Lats: Latitude coordinates of the points of interest (City, Landfills...)
            #Long: Longitude coordinates of the points of interest 
            #region area: Dictionary which defines the limit of the area (cornes latidue and longitudes)
            #patches: The pixels of the overpass
            #coords_plume: Coordinates plume pixels in S5P 'units' calculated before.
        
        lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
        lon_min, lon_max = region_area['lon_min'], region_area['lon_max']
               
        [lon_new, lat_new, index_landfill, index_nearest_landfill,index_nearest_landfill2]=BackSubstractionMethods.rotation_axis(lat, lon, Lats, Long)
    
        #cuadrantes
        area0,area1,area2,area3=[],[],[],[]
        for i in range(len(lat_new)):
            #def cuadrante1
            if lat_new[i]<lat_max and lat_new[i]>lat_new[index_landfill] and lon_new[i]<lon_max and lon_new[i]>lon_new[index_landfill]:
                area0.append(i)
            if lat_new[i]>lat_min and lat_new[i]<lat_new[index_landfill] and lon_new[i]<lon_max and lon_new[i]>lon_new[index_landfill]:
                area1.append(i)
            if lat_new[i]>lat_min and lat_new[i]<lat_new[index_landfill] and lon_new[i]>lon_min and lon_new[i]<lon_new[index_landfill]:
                area2.append(i)
            if lat_new[i]<lat_max and lat_new[i]>lat_new[index_landfill] and lon_new[i]>lon_min and lon_new[i]<lon_new[index_landfill]:
                area3.append(i)    
                
        #only save the pixel near to the landfill
        d=0.5 #distance away (degrees units)
        #check the plume is surrounding the landfill 
        #and delete de plume pixels selected which are far away (d in deg) from the landfill
        areas=[area0,area1,area2,area3]
        area0_new,area1_new,area2_new,area3_new=[],[],[],[]
        areas_new=[area0_new,area1_new,area2_new,area3_new]
        for j,area in enumerate(areas):
            for pixel_idx in area:
                delta_long=lon[pixel_idx]-Long[0]
                delta_lats=lat[pixel_idx]-Lats[0]
                mod=np.sqrt(delta_long**2+delta_lats**2) #module distance
                if mod<=d:
                    areas_new[j].append(pixel_idx)
        
        #calculate the number of plume pixels located in each area
        intersection_pixels=np.zeros(len(areas_new))
        for i, area in enumerate(areas_new):
            set1=set(coords_plume)
            set2=set(area)
            common_pixels=set1.intersection(set2)
            intersection_pixels[i]=len(common_pixels)
        
        #the area selected from the background will be the opposite area to the area which has the max number of plume pixels
        n=np.argmax(intersection_pixels)
        idx_area_selected=n-2
        
        #calculate bck from selected area
        bck=np.mean(var[areas_new[idx_area_selected]])
        var_bck_corrected_wind=var-bck
        
        """
        fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[8,10]  )
        ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
        p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
        p.set_array(var)
        p.set_cmap('coolwarm')
       # p.set_clim(vmin=-30, vmax=30)
        ax.add_collection(p)
        
        
        ax.add_feature(cf.STATES, linewidth=0.2, edgecolor='grey',zorder=3)
        plt.colorbar(p,ax=ax,label='test',orientation='horizontal',shrink=0.6,pad=0.04)

        ax.plot(lon[index_nearest_landfill],lat[index_nearest_landfill],'r*', markersize=5)
        ax.plot(lon[index_landfill],lat[index_landfill],'y*', markersize=5)
        ax.plot(lon[areas_new[idx_area_selected]],lat[areas_new[idx_area_selected]],'kx', markersize=5)
        ax.plot(lon[coords_plume],lat[coords_plume],'bx', markersize=5)
                  
              
        plt.show()
        """
        
        return var_bck_corrected_wind
    