#IMPORTS
#General Modules
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import os

#for plot methane background trends
import datetime 
from brokenaxes import brokenaxes
import NoaaData as ND

#Plot geo maps
import cartopy.crs as ccrs
import cartopy.feature as cf

#import scripts
#ERA5
import VerticalProfileWinds as VPW #for extract the wind velocity and direction from the profile
import SurfacePressure_10_100winds as SPWF #to extract the surface pressure and the wind fields at 10 and 100 meters
#Tropomi
import TropomiData as TD #to extract S5P data
import TropomiPlot as TP #to plot the overpasses with some characteristics, vector fields, colorbars, emission rates
import TropomiFootprint as TF #to calculate the area for each plume pixel
from BckSubstrMethods import BackSubstractionMethods as BSM #different methods to substract the background, striping, down-wind
import TropomiFindPlumePixels as TFPP #to find the plume pixels via thresshold approach

#Plot satellite style TEST
import PlotSatelliteStyle as PSS

#MAIN SCRIPT
#here the emission rate of the overpass is estimated with the data extracted from ERA5 and S5P
#different methdos are applied to substract the background from the overpass

#Select region (Pinto-Madrid-Toledo) corners area of the soundings
region_area= {'lat_min':39.5,'lat_max':41.5, 'lon_min':-5,'lon_max':-2.75 }
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

#variables for the last plot, methane background throughtout the years
years=[]
all_backgrounds=[]
all_date=[]
all_std=[]
length_year=[]
Q_emission_rates=[]
Delta_Q_errors=[]

#for folium satellite plot:
coords_plume_folium=[]
long_bounds_folium=[]
lat_bounds_folium=[]
var_back_corrected_folium=[]
date_folium=[]

#important days 2019-2024, carefully selected
#PBL fixed at 1500m
files={
    2024: [
    ('S5P_OFFL_L2__CH4____20240601T121147_20240601T135317_34378_03_020600_20240607T090505' , 2.1),
    ('S5P_OFFL_L2__CH4____20240616T123110_20240616T141239_34591_03_020600_20240618T045433' , 2.1),
    ('S5P_OFFL_L2__CH4____20240624T114030_20240624T132200_34704_03_020600_20240626T040911' , 2.1),
    ('S5P_OFFL_L2__CH4____20240703T121208_20240703T135338_34832_03_020600_20240705T042848' , 2.5),
    ('S5P_OFFL_L2__CH4____20240704T115305_20240704T133435_34846_03_020600_20240706T040506' , 2.1), #?? Errors?
    ('S5P_OFFL_L2__CH4____20240712T124336_20240712T142506_34960_03_020600_20240714T095018' , 2.3),
    ('S5P_OFFL_L2__CH4____20240723T123644_20240723T141814_35116_03_020600_20240725T064738' , 2.3),
    ('S5P_OFFL_L2__CH4____20240724T121739_20240724T135908_35130_03_020600_20240726T072844' , 2.5),
    ('S5P_OFFL_L2__CH4____20240725T115833_20240725T134003_35144_03_020600_20240728T113012' , 2.0),
    ('S5P_OFFL_L2__CH4____20240816T114403_20240816T132533_35456_03_020600_20240818T040118' , 2.2),
    ('S5P_OFFL_L2__CH4____20240807T125441_20240807T143610_35329_03_020600_20240812T063619' , 2.5)
    ]
    # 2023: [
    # ('S5P_OFFL_L2__CH4____20230624T124634_20230624T142803_29512_03_020500_20230626T050140' , 2.8),
    # ('S5P_OFFL_L2__CH4____20230626T120826_20230626T134956_29540_03_020500_20230628T042345' , 2.3),
    # ('S5P_OFFL_L2__CH4____20230701T121436_20230701T135605_29611_03_020500_20230703T043500' , 2.0),
    # ('S5P_OFFL_L2__CH4____20230713T114834_20230713T133003_29781_03_020500_20230715T040824' , 2.8),
    # ('S5P_OFFL_L2__CH4____20230726T124438_20230726T142608_29966_03_020500_20230729T084827' , 2.5),
    # ('S5P_OFFL_L2__CH4____20230801T123123_20230801T141252_30051_03_020500_20230803T045050' , 2.5),
    # ('S5P_OFFL_L2__CH4____20230803T115307_20230803T133437_30079_03_020500_20230805T042722' , 2.5),
    # ('S5P_OFFL_L2__CH4____20230807T121805_20230807T135934_30136_03_020500_20230911T190701' , 2.5),#?? Errors?
    # ('S5P_OFFL_L2__CH4____20230812T122352_20230812T140521_30207_03_020500_20230814T044517' , 2.2),
    # ('S5P_OFFL_L2__CH4____20230819T115119_20230819T133248_30306_03_020500_20230821T040359' , 2.2)
    # ],
    # 2022: [
    # ('S5P_OFFL_L2__CH4____20220628T121206_20220628T135336_24390_02_020301_20220630T041223' , 2.8),
    # ('S5P_OFFL_L2__CH4____20220708T122437_20220708T140606_24532_02_020301_20220710T043551' , 2.5),
    # ('S5P_OFFL_L2__CH4____20220709T120533_20220709T134703_24546_02_020301_20220711T041613' , 2.1),
    # ('S5P_OFFL_L2__CH4____20220714T121144_20220714T135313_24617_02_020301_20220716T043940' , 2.5),
    # ('S5P_OFFL_L2__CH4____20220726T114547_20220726T132716_24787_03_020400_20220728T040908' , 2.0),
    # ('S5P_OFFL_L2__CH4____20220821T115632_20220821T133801_25156_03_020400_20220825T111049' , 2.3),
    # ('S5P_OFFL_L2__CH4____20220826T120225_20220826T134354_25227_03_020400_20220829T054021' , 1.2) #??
    # ],
    # 2021: [
    # ('S5P_OFFL_L2__CH4____20210629T123700_20210629T141830_19226_01_010400_20210701T054116' , 2.5),
    # ('S5P_OFFL_L2__CH4____20210716T121729_20210716T135858_19467_02_020200_20210718T053735' , 2.8),
    # ('S5P_OFFL_L2__CH4____20210717T115824_20210717T133954_19481_02_020200_20210719T052048' , 2.8), #?? Errors?
    # ('S5P_OFFL_L2__CH4____20210827T122710_20210827T140840_20063_02_020200_20210829T061124' , 2.9),
    # ],
    # 2020: [
    # ('S5P_OFFL_L2__CH4____20200804T120348_20200804T134517_14558_01_010302_20200809T082112' , 2.5),
    # ('S5P_OFFL_L2__CH4____20200826T114902_20200826T133031_14870_01_010302_20200828T053834' , 2.5)
    # ], 
    # 2019: [
    # ('S5P_OFFL_L2__CH4____20190802T120357_20190802T134527_09337_01_010302_20190808T135252' , 2.5)
    # ]
}    

for year, array in files.items():
    Folder_name= 'C:/M-MRS/SP5_Data/'+str(year)+'_CH4'  
    length_year.append(len(os.listdir(Folder_name)))
    years.append(year)
    
    for tupla in array:
        File_name=tupla[0]
        #TROPOMI FILE
        #Load an L2 file (each L2 file contains the measurements taken during one individual orbit).
        F_upwind=2.0
        F_stripe=tupla[1]
        F_clas=2.0
        
        #folder name and file dataset Tropomi
        File_name=str(File_name)
        File = Dataset(f'{Folder_name}/{File_name}/{File_name}.nc')

        #path where save the figs and resutls
        date=File_name[20:24]+'_'+File_name[24:26]+'_'+File_name[26:28]
        year=date[:4]
        month=date[5:7]
        day=date[8:10]
        folder_results="C:/M-MRS/Results/"+year+'/'+date
        if not os.path.exists(folder_results):
            os.makedirs(folder_results)
              
        ############################
        ######      ERA5      ######
        ############################
        
        print('\nERA5')
        
        #WIND FIELDS PROFILE DATA 
        Folder_name_profile='C:/M-MRS/ERA5_Data/Wind_profile/'+year+'/'
        File_name='Pressure_Levels_'+date
        File_name_profile=Folder_name_profile+File_name
                
        FileProfile = Dataset(f'{File_name_profile}.nc')
                
        #SURFACE PRESSURE, 10M AND 100M WIND FIELDS DATA
        Folder_name_wind100='C:/M-MRS/ERA5_Data/Surface_pressure_and_10_100_m_wind/'+year+'/'
        File_name = 'Single_Level_'+date
        File_name=Folder_name_wind100+File_name
        
        File_Sing_Level = Dataset(f'{File_name}.nc')
           
        #function
        figname1=f"Vertical_Wind_Profile"
        figname2="Temp_Humid_O3_Profile"
        path_save1=os.path.join(folder_results,figname1)
        path_save2=os.path.join(folder_results,figname2)
        [wind_u_profile, wind_v_profile, wind_mean, PBL_altitude, delta_U]=VPW.vertical_profiel_wind_field(FileProfile,Lats,Long, path_save1=path_save1,path_save2=path_save2)
                
        #function
        [Pressure, winds_array, lat, lon]=SPWF.surface_pressure_wind_field(File_Sing_Level,Lats,Long,region_area)
        
        # Calculations
        wind_final_data={
            'Wind Profile': (wind_mean, wind_u_profile, wind_v_profile)} 
        
        ##############################
        ######      TROPOMI     ######
        ##############################
        
        print('\nTROPOMI S5P')
          
        #function
        [var,lon,lat,long_bounds,lat_bounds, patches,date]=TD.tropomi_outputs(File, Lats, Long, region_area)
        
        #Plot all togheter
        title=f'S5P - XCH$_4$ measurements and Winds vectors'
        colorbar_label='XCH$_4$ [ppb]'
        figname=f"XCH4_Overpass_Winds_vectors"
        path_save=os.path.join(folder_results,figname)
        graph=TP.tropomi_plot(var, Lats, Long, patches, region_area, title, colorbar_label,
                              quivers_dictoniary=wind_final_data, 
                              path_save=path_save)
              
        #mean value
        mean_density=np.mean(var)
        print('Mean value image:',mean_density,'ppb\n')
        
        #################################
        # reduce background CLASSICAL METHOD
        
        #Find plume pixels

        coords_plume_clas , bckg_mean_clas = TFPP.find_plume_pixels(var,lon,lat,Long,Lats, long_bounds, lat_bounds,F=F_clas)

        #bckg_mean is the background without plume pixels     

        #CALCULATIONS
        var_back_corrected=np.zeros_like(var)
        var_back_corrected=var-bckg_mean_clas
        
        
        #save all backgrounds to study progresion throughout the years
        coords_all_var=np.linspace(0,len(var)-1,len(var)-1)
        coords_no_plume=[int(num) for num in coords_all_var if num not in coords_plume_clas]
        std=np.std(var[coords_no_plume])
        
        
        all_backgrounds.append(bckg_mean_clas)
        all_date.append(date[:10])
        all_std.append(std)
        
        #################################
        # reduce background DESTRIPING METHOD
        # print('destriping')

        bck_sub_strip=BSM.striping(var, lat, lon, Lats, Long, region_area ,patches,long_bounds, lat_bounds,f=F_stripe)
        
        #find plume pixels for image reduced striping method

        coords_plume_striping , bckg_mean_striping = TFPP.find_plume_pixels(bck_sub_strip,lon,lat,Long,Lats, long_bounds, lat_bounds, F=F_stripe)
        
        if coords_plume_striping:
            Percent_plume = abs(np.mean(var[coords_plume_striping])-np.mean(var))/np.mean(var) *100
        else:
            Percent_plume = 0
        
        ################ plot satellite map (folium packg)
        coords_plume_folium.append(coords_plume_striping)
        long_bounds_folium.append(long_bounds)
        lat_bounds_folium.append(lat_bounds)
        var_back_corrected_folium.append(bck_sub_strip)
        date_folium.append(date)
        ################
        
        #################################
        # reduce background UP-WIND METHOD
        # print('up-wind')
        
        # this method takes the plume pixels from the classical method ro estimate where is the location of the plume
        bck_sub_wind=BSM.up_wind(var, lat, lon, Lats, Long, region_area ,patches, coords_plume_clas)
        
        #find plume pixels for image reduced striping method
        
        coords_plume_down_wind , bckg_mean_wind = TFPP.find_plume_pixels(bck_sub_wind,lon,lat,Long,Lats,long_bounds, lat_bounds, F=F_upwind)
        
        ##################################
        
        print('DXCH4 ppb striping: mean:',np.mean(bck_sub_strip[coords_plume_striping]),
                                  'total:',np.sum(bck_sub_strip[coords_plume_striping]))
        print('DXCH4 ppb classical: mean:',np.mean(var_back_corrected[coords_plume_clas]),
                                  'total:',np.sum(var_back_corrected[coords_plume_clas]))
        print('DXCH4 ppb down-wind: mean:',np.mean(bck_sub_wind[coords_plume_down_wind]),
                                  'total:',np.sum(bck_sub_wind[coords_plume_down_wind]))
        
        
        #ESTIMATION OF EMISSION RATES
        #constants 
        molar_mass=16.043 #g/mol methane
        N_A=6.022e23 #mol^-1
        M_d=28.949 #g/mol molecular mass of dry air
        g=9.81 #m/s^2
        V_d=Pressure*N_A/(g*M_d)  #column dry air molecules/m^2
        
        #plume coordinates dictionary for each reduction method
        coords_plume_dictionary={'Classical': [coords_plume_clas,var_back_corrected,F_clas],
                                 'Destriping': [coords_plume_striping,bck_sub_strip,F_stripe],
                                 'Up-Wind': [coords_plume_down_wind,bck_sub_wind,F_upwind]}
        for key, (coord_plume,back_corrected,F) in coords_plume_dictionary.items():
            print('\n',key)
            #pixel area, typically 7x5.5 km at nadir
            area, std_area=TF.footprint_area_computed(long_bounds, lat_bounds, coord_plume)
            
            #sum over the plume pixels  
            ppb=np.array(back_corrected[coord_plume]) #parts per billion (ppb) molecules of ch4 in 1 billion of dry air
            ME=ppb*1e-9*V_d*molar_mass/N_A  #g/m^2
            IME=np.sum(ME*area)   #IME is calculated using eq (7) from amt-11-5673-2018 (g)
            std_ME=np.std(ME)
            delta_IME=area*std_ME + np.sum(ME)*std_area
            print('Sum over the plume pixels (IME): ',IME,'\n')
            
            #typical length scale is calculated using the eq(11) from amt-11-5673-2018
            N=len(coord_plume)
            L=np.sqrt(N*area)  #m longitude units
            
            #Effective wind calculate using eq (12) and fig (4) from amt-11-5673-2018
            #for all wind vectors
            Wind_eff={}
            Q={}
            for key2, (mean,u,v) in wind_final_data.items():
                Wind_eff[key2]=mean
                Q[key2]=Wind_eff[key2]/(L)*IME*3.6*24*365/1000 #kt/yr

                print(key2)
                print('Emission rate (t/day)= ',Q[key2]*1000/365 ,
                      '\nEmission rate (kt/yr)= ',Q[key2],'\n') 
                
            delta_U, delta_IME, delta_L = delta_U , delta_IME , 1/2*np.sqrt(N/area)*std_area #units [m/s], [g], [1/m]
            #delta_U is calculated as the standar deviation for each component (u,v) in the PBL and propagation the error to the mean in the PBL
                #the U_eff is computed using the time, spacial and altitude average. In the last one is calculated the std for the errors.
            #delta_IME is calculated as the error propagation of f=ME*are, delta_ME=std_ME and delta_area=std_area
            #delta_L is calculated with the function L=sqrt(N*area) where delta_N=0 and delta_area=std_area
            
            print(delta_U,delta_IME, delta_L)
            Propagation_errors= abs(IME/L)*delta_U + abs(Wind_eff['Wind Profile']/L)*delta_IME + abs(Wind_eff['Wind Profile']*IME/(L**2))*delta_L  #[g/s]
            C=3.6*24*365/1000
            print('TEST Error Wind (kt/yr):',abs(IME/L)*delta_U*C)
            print('TEST Error IME (kt/yr):',abs(Wind_eff['Wind Profile']/L)*delta_IME*C)
            print('TEST Error Plume length (kt/yr):',abs(Wind_eff['Wind Profile']*IME/(L**2))*delta_L*C)
            print('TEST Total Error (kt/yr):',Propagation_errors*C)
            
            Percent_error = abs(1-(Q['Wind Profile']-Propagation_errors*C)/(Q['Wind Profile']))*100
           
            #plots with Emission rates and save the figures
            text1=['Emission Rate (kt/yr): '+str(round(Q['Wind Profile'],3))+f' ± {round(Propagation_errors*C,3)} (Error %: {round(Percent_error,2)})'+f' (PBL: {PBL_altitude}m)\n'+
                  f'STD Thresshold: {F}\n'+
                  f'Percent Plume-Backg (%): {round(Percent_plume,2)}\n'+
                  f'Mean Raw Value: {round(mean_density,2)} (ppb)']
            
            var2=[lon[coord_plume],lat[coord_plume]] 
            title=r"S5P - XCH$_4$ - XCH$_{4,bck}$"f": {key} method ({date})"
            figname=f"XCH4_{key}_method"
            path_save=os.path.join(folder_results,figname)
            
            if key=='Destriping':
                Q_emission_rates.append(round(Q['Wind Profile'],3))
                Delta_Q_errors.append(round(Propagation_errors*C,3))
               
            colorbar_label=r'$\Delta$XCH$_4$ [ppb]'
            
            graph=TP.tropomi_plot(back_corrected, Lats, Long, patches,
                                  region_area, title, colorbar_label, 
                                  quivers_dictoniary=wind_final_data,
                                  clim=True, var2=var2, path_save=path_save, text1=text1)#, text2=text2)

all_date = [datetime.datetime.strptime(all_date, '%Y-%m-%d') for all_date in all_date]
#plots backgrounds
fig = plt.figure(figsize=(20, 10))

months=['6','7','8']
years=years[::-1]
noaa_date, noaa_bck = ND.noaa_data(years, months)

bax = brokenaxes(
    xlims=(
        (
            datetime.datetime(2019, 6, 1),
            datetime.datetime(2019, 9, 1),
        ),
        (
            datetime.datetime(2020, 6, 1),
            datetime.datetime(2020, 9, 1),
        ),
        (
            datetime.datetime(2021, 6, 1),
            datetime.datetime(2021, 9, 1),
        )
        ,
        (
            datetime.datetime(2022, 6, 1),
            datetime.datetime(2022, 9, 1),
        )
        ,
        (
            datetime.datetime(2023, 6, 1),
            datetime.datetime(2023, 9, 1),
        ),
        (
            datetime.datetime(2024, 6, 1),
            datetime.datetime(2024, 9, 1),
        )
    )
)

bax.plot(noaa_date,noaa_bck,marker='o',color='r',label='NOAA: data month average June-August 2019-2023')
bax.errorbar(all_date,all_backgrounds, yerr=all_std, fmt='o', capsize=5, label='SP5: Overpasses selected June-August 2019-2024')

fig.autofmt_xdate()
[x.remove() for x in bax.diag_handles]
bax.draw_diags()

bax.set_xlabel('Overpass Date',labelpad=60, fontsize=15)
bax.set_ylabel('Background XCH$_4$ [ppb]',labelpad=60, fontsize=15)
bax.legend(loc=2, fontsize=15)
bax.tick_params(axis='x', labelsize=12) 
bax.tick_params(axis='y', labelsize=15)
bax.grid()

#background color
bax.axvspan(datetime.datetime(2024, 6, 27, 0, 0), datetime.datetime(2024, 7, 11, 0, 0), color='red', alpha=0.4)

bax.set_title('Background XCH$_4$ [ppb] NOAA vs TROPOMI: June-August (2019-2024)', pad=20, fontsize=15)
plt.show()

#plots emission rates for destriping method:
fig = plt.figure(figsize=(20, 10))

bax = brokenaxes(
    xlims=(
        # (
        #     datetime.datetime(2019, 6, 1),
        #     datetime.datetime(2019, 9, 1),
        # ),
        # (
        #     datetime.datetime(2020, 6, 1),
        #     datetime.datetime(2020, 9, 1),
        # ),
        # (
        #     datetime.datetime(2021, 6, 1),
        #     datetime.datetime(2021, 9, 1),
        # )
        # ,
        (
            datetime.datetime(2022, 6, 1),
            datetime.datetime(2022, 9, 1),
        )
        ,
        # (
        #     datetime.datetime(2023, 6, 1),
        #     datetime.datetime(2023, 9, 1),
        # ),
        # (
        #     datetime.datetime(2024, 6, 1),
        #     datetime.datetime(2024, 9, 1),
        # )
    )
)

prtr_values_pinto = [12, 12.9, 13.6 ,0.321]
prtr_values_valde = [8.2, 2.5, 2.8 , 1.5]

Krautwurst_value = 13.3 * 24*365/1000 #kt/yr
Krautwurst_value_error=Krautwurst_value*0.26

tu_value=63.428
error_tu_value=5.485

colors = plt.cm.tab10(np.linspace(0, 1, len(prtr_values_pinto)))
for idx, value in enumerate(prtr_values_pinto):
    bax.plot([datetime.datetime(2019+idx, 5, 27, 0, 0),datetime.datetime(2019+idx, 9, 3, 0, 0)],
             [value+prtr_values_valde[idx],value+prtr_values_valde[idx]], color=colors[idx], linestyle='--', label=f'PRTR, Madrid (Pinto + Vald.) {2019+idx}')
 
bax.axhline(y=tu_value, color='purple', linestyle='--', label='Tu et al. (2022): 2018 Sept-Oct campaign')
bax.axhspan(tu_value - error_tu_value, tu_value + error_tu_value, color='purple', alpha=0.4)

bax.axhline(y=Krautwurst_value, color='g', linestyle='--', label='Krautwurst et al. (2024): 2022 August campaign')
bax.axhspan(Krautwurst_value - Krautwurst_value_error, Krautwurst_value + Krautwurst_value_error, color='green', alpha=0.4)

bax.errorbar(all_date,Q_emission_rates, yerr=Delta_Q_errors, fmt='o', capsize=5, label='SP5 data')

fig.autofmt_xdate()
[x.remove() for x in bax.diag_handles]
bax.draw_diags()

bax.set_ylabel('Emission Rates (kt/yr)', labelpad=60, fontsize=15)
bax.set_xlabel('Overpass days', labelpad=60, fontsize=15)
bax.grid()

bax.tick_params(axis='x', labelsize=12) 
bax.tick_params(axis='y', labelsize=15)

#background color
bax.axvspan(datetime.datetime(2024, 6, 27, 0, 0), datetime.datetime(2024, 7, 11, 0, 0), color='red', alpha=0.4)

bax.set_ylim(0,600)

bax.legend(loc=2,fontsize=15)
bax.set_title('Emission Rates [kt/year]: PRTR Spain vs Tu et al. (2022) vs INTA TROPOMI : June-August (2019-2024)',fontsize=15, pad=20)

plt.show()  

########## plt folium map
PSS.plot_satellite_style(coords_plume_folium, long_bounds_folium, lat_bounds_folium, var_back_corrected_folium,date_folium)
        
print('\nEND')

