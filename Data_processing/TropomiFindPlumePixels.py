#Thresshold approach for identificate the plume pixels

#IMPORTS
#General Modules
import numpy as np
import TropomiPlot as TP #to plot the overpasses with some characteristics, vector fields, colorbars, emission rates


f=2 #'detection' limit
def find_plume_pixels(var,lon,lat,Long,Lats, long_bounds, lat_bounds,F=f):
    #Inputs
        #Lats: Latitude coordinates of the points of interest (City, Landfills...)
        #Long: Longitude coordinates of the points of interest
        #F: The 'detection limit' is the number of standard deviations that a pixel has to be above for be plume 

    #Plume detection threshold approach, iteration process to select the plume
   # n=3 #number of iterations
    bckg=np.array(var)
    plume=[]
    #for i in range(n):
    key=True
    n=0
    while key == True:
        n+=1
        #print(n)
        
        bckg_mean=np.mean(bckg)  #background from the total image
        std_bckg=np.std(bckg)
        T=bckg_mean+F*std_bckg  #Thresshold
        plume_new= [x for x in bckg if x > T]
        plume=plume+plume_new
        bckg= [x for x in bckg if x <= T]
        if not plume_new:
            key = False
            
    
    #Plume Coordinates in Tropomi matrix units
    coords_plume=[]
    for i in range(len(plume)):
        coords_plume.append(np.argmin(abs(var-plume[i])))  
        
    d=0.5 #distance away (degrees units)
    #check the plume is surrounding the landfill 
    #and delete de plume pixels selected which are far away (d in deg) from the landfill
    coords_plume_new=[]
    for i in range(len(coords_plume)):
        delta_long=lon[coords_plume[i]]-Long[0]
        delta_lats=lat[coords_plume[i]]-Lats[0]
        mod=np.sqrt(delta_long**2+delta_lats**2) #module distance
        if mod<=d:
            coords_plume_new.append(coords_plume[i])
                 
    #filter to avoid margin plume pixels selection
    """
    if len(coords_plume_new) > 1:
        pixel1, pixel2 =  coords_plume_new[0],  coords_plume_new[0]+1 #take one plume pixel and the following pixel (+1)
        #limit:
        mod = np.sqrt( (lon[pixel1]-lon[pixel2])**2 + (lat[pixel1]-lat[pixel2])**2) #all plume pixels must have a distance below to this to each others
        
        coords_plume_old=coords_plume_new.copy() #to create a different variable
        
        for plume_pixels_i in coords_plume_new:
            check=[]
            
            for plume_pixels_j in coords_plume_old:
                
                if plume_pixels_i ==  plume_pixels_j:
                    continue
                else:
                    mod_plume = np.sqrt( (lon[plume_pixels_i]-lon[plume_pixels_j])**2 + (lat[plume_pixels_i]-lat[plume_pixels_j])**2)
                    
                    if mod_plume >= mod:
                        check.append(True) 
                    else:
                        check.append(False)
            if all(check) == True: 
                if len(coords_plume_old) == 1:
                    continue
                
                if len(coords_plume_old) == 2: #if there are only 2 separeted pixels, then remove only the farthest pixel
                    mod_i = np.sqrt( (lon[plume_pixels_i]-Long[0])**2 + (lat[plume_pixels_i]-Lats[0])**2)
                    mod_j = np.sqrt( (lon[plume_pixels_j]-Long[0])**2 + (lat[plume_pixels_j]-Lats[0])**2)
                    if mod_i<mod_j:
                        coords_plume_old.remove(plume_pixels_j)
                        
                    else:
                        coords_plume_old.remove(plume_pixels_i)
                else:
                    coords_plume_old.remove(plume_pixels_i)
    
        coords_plume_new = coords_plume_old
    """   
        
    #filter to select only plume pixels which must be in contact
    # print(coords_plume_new)
    all_pixels_in_contact=[]
    if len(coords_plume_new) > 2:
        coords_plume_old=coords_plume_new.copy()
        
        for i, pixel_i in enumerate(coords_plume_new):
            corners_pixel_i = set([(round(long_bounds[0][pixel_i],3),round(lat_bounds[0][pixel_i],3)),(round(long_bounds[1][pixel_i],3),round(lat_bounds[1][pixel_i],3)),
                                   (round(long_bounds[2][pixel_i],3),round(lat_bounds[2][pixel_i],3)),(round(long_bounds[3][pixel_i],3),round(lat_bounds[3][pixel_i],3))])
            check=[]
            pixels_contact=[pixel_i]
            for j, pixel_j in enumerate(coords_plume_old):
                if pixel_i==pixel_j:
                    continue
                else:
                    corners_pixel_j = set([(round(long_bounds[0][pixel_j],3),round(lat_bounds[0][pixel_j],3)),(round(long_bounds[1][pixel_j],3),round(lat_bounds[1][pixel_j],3)),
                                           (round(long_bounds[2][pixel_j],3),round(lat_bounds[2][pixel_j],3)),(round(long_bounds[3][pixel_j],3),round(lat_bounds[3][pixel_j],3))])
                    
                    if corners_pixel_i & corners_pixel_j:
                        pixels_contact.append(pixel_j)
                    else:
                        check.append(1)  
            all_pixels_in_contact.append(pixels_contact)
            if np.sum(check)==len(coords_plume_old)-1:
                coords_plume_old.remove(pixel_i)
                   
        coords_plume_new = coords_plume_old   
        
    #filter to avoid individual groups of pixels
    groups=all_pixels_in_contact
    groups_new=[]
    while len(groups_new) != len(groups):
        groups_new = groups.copy()
        for array_i in groups:
            n=0
            for array_j in groups:
                if array_i==array_j:
                    continue
                elif set(array_i) & set(array_j):
                    groups.remove(array_j)
                    groups.append(list(set(array_i).union(set(array_j))))
                    n+=1     
            if n != 0:
                groups.remove(array_i)
            elif n==0 and len(array_i) == 1:
                groups.remove(array_i)
    
    array_not_repeated = set(tuple(arr) for arr in groups)
    groups1=[list(arr) for arr in array_not_repeated]

    len_groups=[]
    for i in range(len(groups1)):
        len_groups.append(len(groups1[i]))
   
    if len_groups:
        len_max=np.max(len_groups)
        idx_max=np.argmin(abs(len_groups-len_max))
        n_pixels_max=len_groups[idx_max]
        
        #we neglect those groups which have fewer pixel than the largest group
        for idx , n_len in enumerate(len_groups):
            if n_len<len_max:
                groups1[idx]=0
            
        
    coords_plume_final=[]
    for group in groups1:
        if group != 0:
            coords_plume_final = coords_plume_final + group
                
    return coords_plume_final , bckg_mean
        

