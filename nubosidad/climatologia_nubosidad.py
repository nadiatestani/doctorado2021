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

#ejemplo:
dset["cldamt_types"].values[0][17] #para cada tipo de nube un array con el porcentaje de ese tipo de nube en cada celda. SOn 18 tipos de nubes
len(dset["cldamt_types"].values[0][7][:,1]) #180 filas
len(dset["cldamt_types"].values[0][7][1,:]) #360 columnas
dset["cldamt_types"].description
"""
<xarray.DataArray 'cloud_type_label' (cloud_type: 18)>
array([b'cumulus_liquid         (680 < PC <= 1025hPa, 0 <= TAU <= 3.55, TC >= 253K)      ',
       b'stratocumulus_liquid   (680 < PC <= 1025hPa, 3.55 < TAU <= 22.63, TC >= 253K)   ',
       b'stratus_liquid         (680 < PC <= 1025hPa, 22.63 < TAU <= 450, TC >= 253K)    ',
       b'cumulus_ice            (680 < PC <= 1025hPa, 0 <= TAU <= 3.55, TC < 253K)       ',
       b'stratocumulus_ice      (680 < PC <= 1025hPa, 3.55 < TAU <= 22.63, TC < 253K)    ',
       b'stratus_ice            (680 < PC <= 1025hPa, 22.63 < TAU <= 450, TC < 253K)     ',
       b'altocumulus_liquid     (440 < PC <= 680hPa, 0 <= TAU <= 3.55, TC >= 253K)       ',
       b'altostratus_liquid     (440 < PC <= 680hPa, 3.55 < TAU <= 22.63, TC >= 253K)    ',
       b'nimbostratus_liquid    (440 < PC <= 680hPa, 22.63 < TAU <= 450, TC >= 253K)     ',
       b'altocumulus_ice        (440 < PC <= 680hPa, 0 <= TAU <= 3.55, TC < 253K)        ',
       b'altostratus_ice        (440 < PC <= 680hPa, 3.55 < TAU <= 22.63, TC < 253K)     ',
       b'nimbostratus_ice       (440 < PC <= 680hPa, 22.63 < TAU <= 450, TC < 253K)      ',
       b'cirrus_liquid          (10 <= PC <= 440hPa, 0 <= TAU <= 3.55, TC >= 253K)       ',
       b'cirrostratus_liquid    (10 <= PC <= 440hPa, 3.55 < TAU <= 22.63, TC >= 253K)    ',
       b'deep_convective_liquid (10 <= PC <= 440hPa, 22.63 < TAU <= 450, TC >= 253K)     ',
       b'cirrus_ice             (10 <= PC <= 440hPa, 0 <= TAU <= 3.55, TC < 253K)        ',
       b'cirrostratus_ice       (10 <= PC <= 440hPa, 3.55 < TAU <= 22.63, TC < 253K)     ',
       b'deep_convective_ice    (10 <= PC <= 440hPa, 22.63 < TAU <= 450, TC < 253K)      '],
      dtype='|S80')
Dimensions without coordinates: cloud_type
Attributes:
    long_name:  Cloud type label
    units:      1
"""


dset["n_orig"] #Mean hourly frequency of occurance of mean cloud amount

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

"""
defino funcion que grafica los campos de una determinada variable en una determinada region
"""

def grafico_campos_nubosidad(paises,provincias,data_list,indice_list,variable,lat_min,lat_max,lon_min,lon_max,unidades_nombre,valor_minimo, valor_maximo, delta_valor,xticks_min,xticks_max, yticks_min, yticks_max,grid,region, ruta_salida, paleta_color):
    """
    Parameters
    ----------
    paises : shapely.geometry.multipolygon.MultiPolygon
        shape con paises a graficar en mapa
    provincias : shapely.geometry.multipolygon.MultiPolygon
        shape con provincias a graficar en mapa
    data_list : list
        lista con data, en cada elemento de la lista hay un NetCDF i.e un xarray.core.dataset.Dataset
    indice_list : float
        indice del elemento de la lista a abrir
    variable : string
        nombre de la variable de los NetCDF a graficar
    lat_min : float
        latitud minima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado
    lat_max : float
        latitud maxima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado
    lon_min : float
        longitud minima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado, usar grados oeste con su signo -
    lon_max : float
        longitud maxima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado, usar grados oeste con su signo -
    unidades_nombre : string
        nombre que se desea para las unidades de la variable, aparece en el grafico
    valor_minimo : float
        valor minimo que toma la escala de la variable
    valor_maximo : float
        valor maximo que toma la escala de la variable
    delta_valor : float
        intervalo entre valores de la escala de la variable
    xticks_min : float
        minimo de longitud que se marca en el grafico en grados oeste con el signo -
    xticks_max : float
        maximo de longitud que se marca en el grafico en grados oeste con el signo -
    yticks_min : float
        minimo de latitud que se marca en el grafico en grados sur con el signo -.
    yticks_max : float
        maximo de latitud que se marca en el grafico en grados sur con el signo -.
    grid : optional
        The default is True.
    region: string
        Nombre de la region
    ruta_salida : str
        Ruta donde se guardan los graficos
    paleta_color: rain (de cero a positivos) / curl (para negativos y positivos)

    Returns
    -------
    Guarda graficos y los muestra.

    """
    #carga librerias necesarias
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cmocean
    
    #limpio graficos
    plt.close()
    
    #selecciona variable
    data=data_list[indice_list]
    variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    variable_data.attrs["units"]=unidades_nombre #cambio el nombre de la unidad

    #selecciono region
    lats=variable_data["lat"][:]
    lons=variable_data["lon"][:]
    lat_lims=[lat_min,lat_max]
    lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
    lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
    lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
    variable_data_subset=variable_data[lat_inds,lon_inds]

    #extraigo mes
    mes=str(data_list[indice_list]["time"].values[0])[5:7]
    #extraigo anio
    anio=str(data_list[indice_list]["time"].values[0])[0:4]
    
    #ploteo
    fig1 = plt.figure(figsize=[12,5],dpi=200)
    ax = fig1.add_subplot(111,projection=ccrs.PlateCarree(central_longitude=0))

    if (paleta_color=="rain"):
        variable_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(valor_minimo, valor_maximo, delta_valor),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': variable_data_subset.units},
                   cmap=cmocean.cm.rain)
    
    if (paleta_color=="curl"):
        variable_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(valor_minimo, valor_maximo, delta_valor),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': variable_data_subset.units},
                   cmap=cmocean.cm.curl_r)

    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.5',linewidth=0.7,alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.4',alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xticklabels(np.arange(xticks_min,xticks_max)[::8])
    plt.xticks(np.arange(xticks_min,xticks_max)[::8])
    ax.set_xlabel("Longitud")

    ax.set_yticklabels(np.arange(yticks_min,yticks_max)[::8])
    plt.yticks(np.arange(yticks_min,yticks_max)[::8])
    ax.set_ylabel("Latitud")

    if (grid==True):
        plt.grid(linestyle="--", alpha=0.3)

    plt.title(variable+" "+region+" "+mes+"/"+anio)
    #plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes+"-"+anio)
    plt.show()


#%%
"""
Grafico los campos de amount de nubosidad (%): cldamt para todos los tiempos y los guardo
"""

import numpy as np

#cargo shape con paises
from shapely.geometry.multipolygon import MultiPolygon
from cartopy.io import shapereader
import geopandas

resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'
shpfilename = shapereader.natural_earth(resolution, category, name) #cargo paises de natural_earth
# cargo el shapefile usando geopandas
df = geopandas.read_file(shpfilename)
# leo los paises que voy a necesitar
paises=MultiPolygon([df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Brazil']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Paraguay']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Uruguay']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Bolivia']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Chile']['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Colombia"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Ecuador"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Venezuela"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Guyana"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Suriname"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Panama"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Costa Rica"]['geometry'].values[0]]) #los paso a multipolygon para poder graficarlos

#cargo shape con provincias de argentina con datos del IGN 
#descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN=geopandas.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias=[None]*24
for i in range(0,24):
    provincias[i]=IGN["geometry"][i]
provincias=MultiPolygon(provincias) #paso a multipolygon para poder ponerlo en mapa


for i in range(0,cantidad_de_datos):
    grafico_campos_nubosidad(paises,provincias,data_list,i,"cldamt",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos","rain")

#grafico_campos_nubosidad(paises,provincias,data_list,5,"cldamt",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos","rain")

#%%
"""

Calculo anomalias climaticas mensuales de variable. Tomo los anios enteros: 1984 a 2016 (32 anios) (DUDA seria hasta 2014?)

1) Calculo media climatica de variable mensual para cada mes, es decir por ejemplo para Enero seria el promedio de la cantidad punto a punto de todos los Eneros del periodo climatologico.

2) A cada mes de cada anio punto a punto le resto la media de ese mes calculada en 1

"""

#%%
"""
Calculo media climatologica por estación y por mes.
Tomo de 1984 a 2016 (32 anios enteros)

Guardo estas cantidades como xarrays y los incorporo a la lista con xarrays por mes.
"""

#hago un array 3d en donde cada capa es un campo mensual y calculo la media
def media_mensual(data_list,variable,mes):
    """

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio
    variable : str
        nombre variable
    mes : str
        numero de mes dos digitos

    Returns
    -------
    Array con media del mes seleccionado

    """
    import numpy as np
    n=1
    arr_3D=np.empty((1,180,360))
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes):
            variable_data=[data_list[i][variable].mean("time", keep_attrs=True).values]
            arr_3D=np.concatenate([arr_3D,variable_data])
            n=n+1
            arr_3D=np.reshape(arr_3D,(n,180,360))
    arr_3D=arr_3D[1:np.shape(arr_3D)[0],:,:]
    media_mensual=np.mean(arr_3D,axis=0)
    return(media_mensual)

#hago la media para cada mes y lo incorporo a una lista como xarray

def media_mensual_xarray(data_list,variable,mes):
    import numpy as np
    import xarray as xr
    data=media_mensual(data_list, variable, mes)
    lats=data_list[2][variable].mean("time", keep_attrs=True)["lat"][:].values
    lons=data_list[2][variable].mean("time", keep_attrs=True)["lon"][:].values
    #dims=["lat", "lon"]
    coords=[("lat", lats),("lon", lons)]
    #coords=data_list[0][variable].coords
    xarray_salida=xr.DataArray(data, coords=coords)
    xarray_salida.name=variable+ " media mensual "+ mes
    xarray_salida.attrs['units']="%"
    xarray_salida.attrs['mes']=mes
    return(xarray_salida)

media_mensual_xarray(data_list,"cldamt","09")

#lista_media_mensual=...
#SEGUIR DESDE ACA
#%%
#hago un array 3d en donde cada capa es un campo mensual y calculo la media
def media_mensual(data_list,variable,mes):
    """

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio
    variable : str
        nombre variable
    mes : str
        numero de mes dos digitos

    Returns
    -------
    Array con media del mes seleccionado

    """
    import numpy as np
    n=1
    arr_3D=np.empty((1,180,360))
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes):
            variable_data=[data_list[i][variable].mean("time", keep_attrs=True).values]
            arr_3D=np.concatenate([arr_3D,variable_data])
            n=n+1
            arr_3D=np.reshape(arr_3D,(n,180,360))
    arr_3D=arr_3D[1:np.shape(arr_3D)[0],:,:]
    media_mensual=np.mean(arr_3D,axis=0)
    return(media_mensual)

#calculo la anomalia mensual restandole a cada mes la media de ese mes y lo agrego al xarray correspondiente de la lista data_list con nombre "media_climatologica_"+variable
def anomalia_mensual(data_list,variable,mes):
    """

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, cargar la lista para modificar
    variable : str
        nombre variable
    mes : str
        numero de mes dos digitos

    Returns
    -------
    la lista con la variable de anomalia agregada al xarray

    """
    import xarray as xr
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes):
            variable_data=[data_list[i][variable].mean("time", keep_attrs=True).values]
            anom=variable_data-media_mensual(data_list,variable,mes)
            anom_dataarray=xr.DataArray(data=anom,dims=["time","lat","lon"])
            data_list[i]=data_list[i].assign(variable_anom=anom_dataarray)
            data_list[i]=data_list[i].rename({"variable_anom":"anomalia_mensual_"+variable})
    return(data_list)

data_list_modificado=data_list.copy()
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
for i in range(0,12):
    data_list_modificado=anomalia_mensual(data_list_modificado,"cldamt",meses[i])


for i in range(0,cantidad_de_datos):
    grafico_campos_nubosidad(paises,provincias,data_list_modificado,i,"anomalia_mensual_cldamt",-39,-16,-64,-31,"%",-60,65,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_anomalias_mensuales_campos","curl")
    
#%%

"""
calculo la media mensual de una determinada variable en determinada region 
grafico una serie temporal de todo el periodo, y por mes
"""
#defino funcion que calcula media mensual de una determinada variable en determinada region

def media_espacial(data_list,indice_list,variable,lat_min,lat_max,lon_min,lon_max):
    """
    Calcula media mensual de una determinada variable en determinada region

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, cargar la lista para modificar
    indice_list: int
        indice de elemento de la lista en el que busca los datos
    variable : str
        nombre variable
    lat_min : float
        latitud minima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado
    lat_max : float
        latitud maxima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado
    lon_min : float
        longitud minima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado, usar grados oeste con su signo -
    lon_max : float
        longitud maxima a seleccionar de las variables, obs: van con decimales 0.5 y cada 1 grado, usar grados oeste con su signo -

    Returns float
    -------
    media espacial de la variable seleccionada en la region seleccionada
    """
    import numpy as np
    
    data=data_list[indice_list]
    variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    
    #selecciono region
    lats=variable_data["lat"][:]
    lons=variable_data["lon"][:]
    lat_lims=[lat_min,lat_max]
    lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
    lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
    lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
    variable_data_subset=variable_data[lat_inds,lon_inds]
    
    #calculo media espacial del subset de la variable data
    media_espacial=np.mean(variable_data_subset.values,axis=(0,1))
    return(media_espacial)

#armo funcion que devuelve data frame donde primera columna sea la fecha y la segunda columna sea la media espacial
def media_espacial_df(data_list,variable,lat_min,lat_max,lon_min,lon_max):
    media_espacial_df=pd.DataFrame(columns=["fecha","Media_espacial_"+variable])
    for i in range(0,len(data_list)):
        #extraigo mes
        mes=str(data_list[i]["time"].values[0])[5:7]
        #extraigo anio
        anio=str(data_list[i]["time"].values[0])[0:4]
        #¢alculo media
        media_espacial_i=media_espacial(data_list,i,variable,lat_min,lat_max,lon_min,lon_max)
        media_espacial_df=media_espacial_df.append({"fecha": mes+"-"+anio, "Media_espacial_"+variable :media_espacial_i},ignore_index=True)
        
    return(media_espacial_df)

cldamt_media_espacial_df=media_espacial_df(data_list,"cldamt",-60,-31,-35,-18)

import matplotlib.pyplot as plt


fig1, ax = plt.subplots(figsize=[12,6],dpi=200)
plt.plot(cldamt_media_espacial_df["fecha"],cldamt_media_espacial_df["Media_espacial_cldamt"])
major_ticksx=np.arange(6,len(cldamt_media_espacial_df["fecha"]),len(cldamt_media_espacial_df["fecha"])/34)
ax.set_xticks(major_ticksx)
ax.grid(alpha=0.3)
ax.tick_params(axis='x',direction='out',bottom=True,labelrotation=45, labelsize=10,pad=1.5)
ax.set_xlabel("Fecha", size=10)
ax.set_ylabel("cldamt media %", size=10)
plt.title("Media espacial mensual cldmt region 1")
plt.savefig("/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos/serie", dpi=140)
