import folium
import branca.colormap as cm
import numpy as np
from folium import GeoJson
from folium.plugins import FloatImage

# Create a map object, centered on a specific location (latitude, longitude)
location = [40.41831, -3.70275]  

Lats=np.array([40.257118, 40.336427, 40.457619, 40.41831])
Long=np.array([-3.637550, -3.590375, -3.363080, -3.70275])
names=['Pinto','Valdemingomez','Alcalá', 'Madrid City']
# Función para agregar una cuadrícula de longitudes y latitudes
def agregar_cuadricula(mapa, lat_start, lat_end, lon_start, lon_end, interval):
    # Agregar líneas de latitud
    for lat in range(int(lat_start), int(lat_end) + 1, interval):
        folium.PolyLine(
            locations=[[lat, lon_start], [lat, lon_end]],
            color='gray',
            weight=1.5,
            opacity=1
        ).add_to(mapa)
    
    # Agregar líneas de longitud
    for lon in range(int(lon_start), int(lon_end) + 1, interval):
        folium.PolyLine(
            locations=[[lat_start, lon], [lat_end, lon]],
            color='gray',
            weight=1.5,
            opacity=1
        ).add_to(mapa)


def plot_satellite_style(plumes_coordinates,long_bounds,lat_bounds,var,date):
    #coordinates of pixels plumes
    #coordiantes of their coorners
    #variables (value)

    my_map = folium.Map(location=location, zoom_start=10)#, tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', attr='Esri')
    
    #add different map stiles
    #1
    folium.TileLayer(
        tiles='CartoDB positron',
        attr='Carto',
        name='Light mode',
        overlay=False,
        control=True).add_to(my_map)  # Estilo claro (calles)
    
    folium.TileLayer(
        tiles='CartoDB dark_matter',
        attr='CartoD',
        name='Dark mode',
        overlay=False,
        control=True).add_to(my_map)  # Estilo oscuro
    #2
    folium.TileLayer(
        tiles='https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}',
        attr='Google',
        name='Google Satellite',
        overlay=False,
        control=True).add_to(my_map)

    #3
    folium.TileLayer(
        tiles='https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}',
        attr='Esri',
        name='Esri Satellite',
        overlay=False,
        control=True).add_to(my_map)
    

    min_value = 10
    max_value = 40
    
    for index , plume_coords in enumerate(plumes_coordinates):
        #poligon data
        polygon_data=[]
        for plume_idx in plume_coords:
            polygon_data.append({"coordinates":[[lat_bounds[index][0][plume_idx],long_bounds[index][0][plume_idx]],
                                              [lat_bounds[index][1][plume_idx],long_bounds[index][1][plume_idx]],
                                              [lat_bounds[index][2][plume_idx],long_bounds[index][2][plume_idx]],
                                              [lat_bounds[index][3][plume_idx],long_bounds[index][3][plume_idx]]
                                              ],
                                "value": var[index][plume_idx]
                                    })
            
        #create a color map
        colormap = cm.LinearColormap(colors=['#ffcccc', '#ff6666', '#ff3333', '#990000'], vmin=min_value, vmax=max_value)
        colormap.caption = "Values ppb"
        
        colormap.width = 600  # Cambia el ancho del colorbar
        colormap.height = 55  # Cambia la altura del colorbar
        
        #layers
        layer = folium.FeatureGroup(name = 'Date: '+str(date[index]),show=False)
        
        #plot the polygons in the layers
        for poly in polygon_data:
            folium.Polygon(
                locations=poly['coordinates'],
                color=colormap(poly['value']),   # Get color based on value
                fill=True,
                fill_opacity=0.7,
                popup=f"Value: {poly['value']}"
            ).add_to(layer)
            
        #add the layer to the main map
        layer.add_to(my_map)
        
   
    #add the color map to the layer
    colormap.add_to(my_map)   
     
    # Add markers to the map
    for idx , name in enumerate(names):
        folium.Marker([Lats[idx], Long[idx]], popup=name).add_to(my_map)
        
    #add provinces limits
    geojson_file = 'georef-spain-provincia.geojson'
    GeoJson(
        geojson_file,
        name="Límites de las provincias",
        style_function=lambda x: {
            'fillColor': 'lightblue',  # Color de relleno
            'color': 'black',            # Color del borde
            'weight': 1.5,                # Grosor del borde
            'fillOpacity': 0.0          # Opacidad del relleno
        }
    ).add_to(my_map)
      
    
    #Agregar cuadrícula al mapa (ChatGPT)
    agregar_cuadricula(my_map, lat_start=35, lat_end=45, lon_start=-9, lon_end=1, interval=1)
    
    # Agrega un control de capas al mapa
    folium.LayerControl().add_to(my_map)
    
    # Save the map as an HTML file
    my_map.save("map.html")
    
