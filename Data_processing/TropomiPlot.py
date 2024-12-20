#TROPOMI SP5 PLOT DATA

#IMPORTS
#General Modules
import numpy as np
import matplotlib.pyplot as plt

#Plot geo maps
import cartopy.crs as ccrs
import cartopy.feature as cf

#For plotting S5P soundings 
from matplotlib.collections import PatchCollection

QUIVERS_DICTIONARY={}
CLIM=False
VAR2=[]
PATH_SAVE=[]
TEXT1=[]
TEXT2=[]

def tropomi_plot(
        var, Lats, Long, patches, region_area, title, colorbar_label,
        quivers_dictoniary=QUIVERS_DICTIONARY, clim=CLIM, var2=VAR2,
        path_save=PATH_SAVE, text1=TEXT1, text2=TEXT2
        ):
    #Inputs
        #Lats: Latitude coordinates of the points of interest (City, Landfills...)
        #Long: Longitude coordinates of the points of interest
        #patches: The pixels of the overpass
        #region area: Dictionary which defines the limit of the area (cornes latidue and longitudes)
        #var2: for plume pixels, introduce a matrix with longs and lats of the pixels
        #quivers_dictionary: a dictionary where the information of the wind direction vectors are saved
        #clim: vmin and vmax for the colormap and colorbar
        #path_save: for save figure in path
        #text1: text dedicated to draw the emission rates on the plot
        #text2: text dedicated to draw the effective wind speed used to calculate de emission rate

    lat_min, lat_max = region_area['lat_min'], region_area['lat_max']
    lon_min, lon_max = region_area['lon_min'], region_area['lon_max']
    
    #Plot features
    fig, ax = plt.subplots(1, 1, subplot_kw=dict(projection=ccrs.PlateCarree()), figsize=[12,14] )
    ax.set_extent([lon_min,lon_max,lat_min,lat_max], crs=ccrs.PlateCarree())
    
    #Plot soundings
    p = PatchCollection(patches,transform=ccrs.PlateCarree(),zorder=0)
    p.set_array(var)
    p.set_cmap('coolwarm')
    ax.add_collection(p)
           
    #Coordinate-Grid
    grid = ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,linewidth=1,linestyle='--', color='gray', alpha=0.5)
    grid.top_labels = False
    grid.right_labels = False
    grid.xlabel_style = {'size': 12}
    grid.ylabel_style = {'size': 12}
    
    #Plot features
    ax.add_feature(cf.COASTLINE,linewidth=1,zorder=4)
    ax.add_feature(cf.BORDERS,linewidth=1,zorder=3)
    ax.add_feature(cf.STATES, linewidth=1, edgecolor='grey',zorder=3)       
    
    #Madrid City and Facilities
    ax.plot(Long[:-1],Lats[:-1],'k*', markersize=13, transform=ccrs.Geodetic(), label='Landfill Facilities')
    ax.plot(Long[-1],Lats[-1],'r*', markersize=13, transform=ccrs.Geodetic(),label='Madrid City')
    
    if quivers_dictoniary:
        colors = plt.cm.summer(np.linspace(0, 1, len(quivers_dictoniary)))
        i=0
        for key, (mean,u,v) in quivers_dictoniary.items():
            plt.quiver(*[Long[0],Lats[0]],u/mean,v/mean,angles='xy', scale_units='xy', scale=5, color=colors[i], label='Wind'+': '+str(round(mean,3))+' m/s') 
            i+=1
    if clim==True:
        p.set_clim(vmin=-30, vmax=30)
    if var2:
        ax.plot(var2[0],var2[1],'bx', markersize=10)
    if text1: #(0.138, 0.645)
        plt.text(0.138, 0.725, text1[0], fontsize=13, color='black', transform=plt.gcf().transFigure, bbox={'facecolor':'white', 'pad':5, 'alpha':0.5})
    if text2:
        plt.text(0.51, 0.79, text2[0], fontsize=13, color='black', transform=plt.gcf().transFigure, bbox={'facecolor':'white', 'pad':5, 'alpha':0.5})
        
    plt.title(title)
    plt.colorbar(p,ax=ax,label=colorbar_label,orientation='horizontal',shrink=0.6,pad=0.04)
    ax.legend(fontsize=13,loc='upper left')
    
    if path_save:
        plt.savefig(path_save)
        
    plt.show()
    plt.close()
