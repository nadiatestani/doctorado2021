#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climatologia de nubosidad con datos de ISCCP-H HGM 

Los datos vienen dados como Promedios mensuales con una resolucion espacial de 1º, estan en formato netCDF CF complaint y se descargaron directamente o con THREDDS. Los detalles de los datos estan en este link: https://www.ncdc.noaa.gov/isccp/isccp-data-access/isccp-basic-data y la metadata aca: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00956 

Los promedios mensuales HGM surgen de hacer el promedio mensual de los datos ISCCP-H HGH que son los datos del promedio cada 3 horas para cada hora sinoptica. 

Los datos empiezan en el mes 07 del año 1983 y se extienden hasta el mes 06 del año 2017

Para descargar la informacion use la opcion de descarga directa con este comando en la terminal (fuente: https://stackoverflow.com/questions/6827459/can-i-use-wget-to-download-multiple-files-from-linux-terminal) :

    >wget -r -l1 -A.nc https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access/isccp-basic/hgm/


* La primera variable con la que trabajo es cldamt: cloud amount (en porcentaje)

Uso este tutorial https://carpentrieslab.github.io/python-aos-lesson/02-visualisation/index.html para construir el codigo
"""

#%% 
"""
Visualizo la informacion:
    Abro uno de los archivos .nc de nubosidad mensual
    Veo las variables
"""

import xarray as xr
nc_ruta="Documentos/Doctorado/datos/nubosidad/ISCCP-H_HGM"
nc_name="ISCCP-Basic.HGM.v01r00.GLOBAL.1986.07.99.9999.GPC.10KM.CS00.EA1.00.nc"
dset=xr.open_dataset(nc_ruta+"/"+nc_name)
print(dset)

"""
Vemos:
    lon (longitudes) va de 0.5 a 359.5 (son 360 longitudes cada 1 grado)
    lat (latitudes) va de -89.5 a 89.5 (son 180 latitudes cada 1 grado)
    time (tiempo) hay un solo tiempo que se corresponde al dia 15 del mes sobre el cual se hace el promedio mensual
    Data variables (algunas):
        
        cldamt: mean cloud amount (en porcentaje)
        cldamt_dist: Average frequency of occurrence of monthly mean cloud amount at each time of day UTC
        pd: Mean cloud pressure
        tc: Mean cloud temperature
        sigma_tc_time: cloud-top temperature (TC) mean standard deviation over time
        cldamt_ir:  Average cloud amount detected by IR threshold regardless of VIS threshold
        cldamt_types: Mean cloud amount for cloud types (%) Cloud detected by either IR or VIS threshold, type determined by cloud top location adjusted for optically thinner clouds and optical thickness for liquid or ice as determined by cloud top temperature

"""
#%%
"""
Abro los datos y los acomodo en una lista. Cada lista se corresponde al dset de cada NetCDF (uno por cada mes por cada anio)
"""
import xarray as xr
import pandas as pd

fecha_inicio="1983-07-15"
fecha_final="2017-07-15"

datelist = pd.date_range(start="1983-07-15",end="2017-07-15",freq="M")
cantidad_de_datos=len(datelist)

#armo lista de nombres
nc_name_list=[None]*cantidad_de_datos
for i in range(0,cantidad_de_datos):
    nc_name_list[i]="ISCCP-Basic.HGM.v01r00.GLOBAL."+str(datelist[i])[0:4]+"."+str(datelist[i])[5:7]+".99.9999.GPC.10KM.CS00.EA1.00.nc"

nc_ruta="Documentos/Doctorado/datos/nubosidad/ISCCP-H_HGM"

data_list=[None]*cantidad_de_datos
for i in range(0,cantidad_de_datos):
    data_list[i]=xr.open_dataset(nc_ruta+"/"+nc_name_list[i])

#%%

#import cartopy.crs as ccrs
#import matplotlib.pyplot as plt
#import numpy as np
#import cartopy.feature as feature
#para graficar cargo libreria de paleta de colores rain
#import cmocean
#para poner pais
#from cartopy.io import shapereader
#import geopandas
#para poner provincia
#from shapely.geometry.multipolygon import MultiPolygon


#veo un primer grafico

cldamt_data=dset["cldamt"].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
cldamt_data.attrs["units"]="%" #cambio el nombre de la unidad

#selecciono region
lats=cldamt_data["lat"][:]
lons=cldamt_data["lon"][:]
lat_lims=[-39,-16] #lean Los datos van de 1 en 1 pero con .5 osea: -38.5 -16.5
lon_lims=[360-64,360-31] #lean 360-64 (64 O) 360-31 (31 O) 
#lat_lims=[-54,-16] #
#lon_lims=[360-80,360-31] #lean 360-64 (64 O) 360-31 (31 O)
lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]

cldamt_data_subset=cldamt_data[lat_inds,lon_inds]


#cargo paises
#descargo mapa de natural_earth
# get country borders
resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'
shpfilename = shapereader.natural_earth(resolution, category, name)
# read the shapefile using geopandas
df = geopandas.read_file(shpfilename)
# read the argentina borders
argentina = df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0]
brasil = df.loc[df['ADMIN'] == 'Brazil']['geometry'].values[0]
paraguay = df.loc[df['ADMIN'] == 'Paraguay']['geometry'].values[0]
uruguay = df.loc[df['ADMIN'] == 'Uruguay']['geometry'].values[0]
bolivia = df.loc[df['ADMIN'] == 'Bolivia']['geometry'].values[0]
chile = df.loc[df['ADMIN'] == 'Chile']['geometry'].values[0]
peru = df.loc[df['ADMIN'] == "Peru"]['geometry'].values[0]

paises=MultiPolygon([argentina,brasil,paraguay,uruguay,bolivia,chile,peru])

#lo hago con datos del IGN 
#descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN=geopandas.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias=[None]*24
for i in range(0,24):
    provincias[i]=IGN["geometry"][i]
provincias=MultiPolygon(provincias) #paso a multipoligonos para poder ponerlo en mapa


#ploteo
fig1 = plt.figure(figsize=[12,5])

ax = fig1.add_subplot(111,projection=ccrs.PlateCarree(central_longitude=0))

cldamt_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(0, 101, 10),
                   extend='max',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': cldamt_data_subset.units},
                   cmap=cmocean.cm.rain)

ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.5',linewidth=0.7,alpha=0.8)

ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.4',alpha=0.8)

ax.coastlines(color='0.3')

ax.set_xticklabels(np.arange(-60,-31)[::4])
plt.xticks(np.arange(-60,-31)[::4])
ax.set_xlabel("Longitud")

ax.set_yticklabels(np.arange(-38,-16)[::4])
plt.yticks(np.arange(-38,-16)[::4])
ax.set_ylabel("Latitud")

plt.grid(linestyle="--", alpha=0.3)

plt.title("Porcentaje de nubosidad anio y mes")
plt.show()


#%%
"""
Abro archivos nc y armo dataframe con todos los datos ?
"""

#%%
"""
Armo funcion que:
    Abre archivos nc
    Selecciona region
    Selecciona variables
"""



