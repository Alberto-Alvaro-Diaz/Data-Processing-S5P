#FOOTPRINT AREA COMPUTE

#IMPORTS
#General Modules
import numpy as np

def footprint_area_computed(long_bounds, lat_bounds, plume_coords):
        #Inputs
            #long_bounds: Longitud cornes coordinates of the each overpass pixel
            #lat_bounds: Latitude cornes coordinates of the each overpass pixel
            #plume_coords: Coordinates plume pixels in S5P 'units' calculated before.
    
    #compute pixel areas. The long and lat coordinates of each corner for each pixel plume are used.
    earth_radio=6378e3 #earth radio in meters
    
    #Other way, using the the corner's coordinate
    lat_x=lat_bounds[0][plume_coords]-lat_bounds[1][plume_coords]  #left bottom corner - right bottom corner
    lon_x=long_bounds[0][plume_coords]-long_bounds[1][plume_coords] 
    mod_x=np.array(np.sqrt(lat_x**2+lon_x**2))*earth_radio*np.pi/180
    
    lat_y=lat_bounds[0][plume_coords]-lat_bounds[3][plume_coords] #left bottom corner - left top corner
    lon_y=long_bounds[0][plume_coords]-long_bounds[3][plume_coords]
    mod_y=np.array(np.sqrt(lat_y**2+lon_y**2))*earth_radio*np.pi/180
    area=mod_x*mod_y
    #print('Area for each plume pixel (m^2): ', area)
    print('Pixel footprint (km x km): ', np.mean(mod_x*1e-3),' x ', np.mean(mod_y*1e-3),'\n')
        
    return np.mean(area), np.std(area)

'''   
    #plot to plume an background selection
    fig, ax = plt.subplots( 1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[8,10]  )
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
    p.set_array(var)
    p.set_cmap('coolwarm')
    ax.add_collection(p)
    ax.add_feature(cf.STATES, linewidth=0.2, edgecolor='grey',zorder=3)
    
    #plot plumes(blue) and bakcground (green) pixels selected
    ax.plot(methane_arr[maxs,1],methane_arr[maxs,2],'b*', markersize=5)
    ax.plot([lon0[maxs], lon1[maxs]],[lat0[maxs],lat1[maxs]], linewidth=1 ,color='blue')
    ax.plot([lon1[maxs], lon2[maxs]],[lat1[maxs],lat2[maxs]], linewidth=1,color='blue')
    ax.plot([lon2[maxs], lon3[maxs]],[lat2[maxs],lat3[maxs]], linewidth=1,color='blue')
    ax.plot([lon3[maxs], lon0[maxs]],[lat3[maxs],lat0[maxs]], linewidth=1,color='blue')

    plt.show()
'''
