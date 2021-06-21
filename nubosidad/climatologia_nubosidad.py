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
Calculos estadisticos básicos
"""
#%%
"""
Calculo media y desvio climatologica por mes.
Tomo de 1984 a 2016 (32 anios enteros)

Guardo estas cantidades como xarrays y los incorporo a la lista con xarrays por mes.
"""

#hago un array 3d en donde cada capa es un campo mensual y calculo la media y desvio, lo acomodo como xarray
def media_desvio_mensual(data_list,variable,mes, unidad):
    """
    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, tomar lista en el rango de anios que se desean
    variable : str
        nombre variable
    mes : str
        numero de mes dos digitos

    Returns
    -------
    [0] xarray con media del mes seleccionado
    [1] xarray con desvio del mes seleccionado

    """
    import numpy as np
    import xarray as xr
    
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
    desvio_mensual=np.std(arr_3D,axis=0)
    
    #cargo lats y lons para armar xarray
    lats=data_list[0][variable].mean("time", keep_attrs=True)["lat"][:].values
    lons=data_list[0][variable].mean("time", keep_attrs=True)["lon"][:].values
    coords=[("lat", lats),("lon", lons)]
    
    #salida media
    xarray_media=xr.DataArray(media_mensual, coords=coords)
    xarray_media.name=variable+ " media mensual "+ mes
    xarray_media.attrs['units']=unidad
    xarray_media.attrs['mes']=mes
    
    #salida desvio
    xarray_desvio=xr.DataArray(desvio_mensual, coords=coords)
    xarray_desvio.name=variable+ " desvio mensual "+ mes
    xarray_desvio.attrs['units']=unidad
    xarray_desvio.attrs['mes']=mes
    
    return([xarray_media,xarray_desvio])

#hago la media y el desvio para cada mes y lo incorporo a dos listas como xarrays, lo hago con anios enteros: 1984-2016
data_list_1984_2016=data_list[6:402]

#armo lista con medias mensuales de cldamt
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
media_mensual_cldamt_list=[None]*12
for i in range(0,12):
    media_mensual_cldamt_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt",meses[i],"%")[0]

#armo lista con desvios mensuales de cldamt
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
desvio_mensual_cldamt_list=[None]*12
for i in range(0,12):
    desvio_mensual_cldamt_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt",meses[i],"%")[1]

#%%
"""
Calculo media climatologica por estación 
Tomo de 1984 a 2016 (32 anios enteros), para que haya estaciones enteras tomo de Diciembre 1983 a Noviembre 2016

Guardo estas cantidades como xarrays y los incorporo a la lista con xarrays por mes.
"""

#hago un array 3d en donde cada capa es un campo mensual y calculo la media
def media_desvio_trimestral(data_list,variable,mes1,mes2,mes3,unidad):
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
    [0] xarray con media del mes seleccionado
    [1] xarray con desvio del mes seleccionado

    """
    import numpy as np
    import xarray as xr
    n1=1
    arr_3D1=np.empty((1,180,360))
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes1):
            variable_data1=[data_list[i][variable].mean("time", keep_attrs=True).values]
            arr_3D1=np.concatenate([arr_3D1,variable_data1])
            n1=n1+1
            arr_3D1=np.reshape(arr_3D1,(n1,180,360))
    arr_3D1=arr_3D1[1:np.shape(arr_3D1)[0],:,:]
    
    n2=1
    arr_3D2=np.empty((1,180,360))
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes2):
            variable_data2=[data_list[i][variable].mean("time", keep_attrs=True).values]
            arr_3D2=np.concatenate([arr_3D2,variable_data2])
            n2=n2+1
            arr_3D2=np.reshape(arr_3D2,(n2,180,360))
    arr_3D2=arr_3D2[1:np.shape(arr_3D2)[0],:,:]
    
    n3=1
    arr_3D3=np.empty((1,180,360))
    for i in range(0,len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7]==mes3):
            variable_data3=[data_list[i][variable].mean("time", keep_attrs=True).values]
            arr_3D3=np.concatenate([arr_3D3,variable_data3])
            n3=n3+1
            arr_3D3=np.reshape(arr_3D3,(n3,180,360))
    arr_3D3=arr_3D3[1:np.shape(arr_3D3)[0],:,:]
    
    arr_3D_3meses=np.concatenate([arr_3D1,arr_3D2,arr_3D3])
    media_trimestral=np.mean(arr_3D_3meses,axis=0)
    desvio_trimestral=np.std(arr_3D_3meses,axis=0)
    
    #cargo lats y lons para armar xarrays
    lats=data_list[0][variable].mean("time", keep_attrs=True)["lat"][:].values
    lons=data_list[0][variable].mean("time", keep_attrs=True)["lon"][:].values
    coords=[("lat", lats),("lon", lons)]
    
    #xarray media
    xarray_media=xr.DataArray(media_trimestral, coords=coords)
    xarray_media.name=variable+ " media trimestral "+ mes1+"-"+mes2+"-"+mes3
    xarray_media.attrs['units']=unidad
    xarray_media.attrs['meses']=[mes1,mes2,mes3]
    
    #xarray media
    xarray_desvio=xr.DataArray(desvio_trimestral, coords=coords)
    xarray_desvio.name=variable+ " desvio trimestral "+ mes1+"-"+mes2+"-"+mes3
    xarray_desvio.attrs['units']=unidad
    xarray_desvio.attrs['meses']=[mes1,mes2,mes3]
    
    return([xarray_media,xarray_desvio])


data_list_12_1983_11_2016=data_list[5:401]

#armo lista con medias trimestrales
media_trimestral_cldamt_list=[None]*4
media_trimestral_cldamt_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","03","04","05","%")[0]
media_trimestral_cldamt_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","06","07","08","%")[0]
media_trimestral_cldamt_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","09","10","11","%")[0]
media_trimestral_cldamt_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","12","01","02","%")[0]

#armo lista con desvios trimestrales
desvio_trimestral_cldamt_list=[None]*4
desvio_trimestral_cldamt_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","03","04","05","%")[1]
desvio_trimestral_cldamt_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","06","07","08","%")[1]
desvio_trimestral_cldamt_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","09","10","11","%")[1]
desvio_trimestral_cldamt_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt","12","01","02","%")[1]

#%%
"""
Defino funcion que grafica climatologia mensual de alguna variable definida previamente en una determinada region
"""

def grafico_campos_climatologia_nubosidad(paises,provincias,data_list_climatologia,indice_list,variable,climatologia_tipo,lat_min,lat_max,lon_min,lon_max,unidades_nombre,valor_minimo, valor_maximo, delta_valor,xticks_min,xticks_max, yticks_min, yticks_max,grid,region, ruta_salida, paleta_color,espacio_entre_lat_lon,orientacion):
    """
    Parameters
    ----------
    paises : shapely.geometry.multipolygon.MultiPolygon
        shape con paises a graficar en mapa
    provincias : shapely.geometry.multipolygon.MultiPolygon
        shape con provincias a graficar en mapa
    data_list_climatologia : list
        lista con data climatologica, en cada elemento de la lista hay un NetCDF i.e un xarray.core.dataset.Dataset
    indice_list : float
        indice del elemento de la lista a abrir
    variable : string
        nombre de la variable de los NetCDF a graficar
    climatologia_tipo: string
        "mensual" o "trimestral" 
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
    paleta_color: rain (de cero a positivos) / curl (para negativos y positivos) /matter (de cero a positivos en rosas)
    espacio_entre_lat_lon: float
        4 para region chica (menos separacion), 8 para region grande (mas separacion)
    orientacion: str
        "H": horizontal "V": vertical
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
    variable_data=data_list_climatologia[indice_list]
    #variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    variable_data.attrs["units"]=unidades_nombre #cambio el nombre de la unidad

    #selecciono region
    lats=variable_data["lat"][:]
    lons=variable_data["lon"][:]
    lat_lims=[lat_min,lat_max]
    lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
    lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
    lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
    variable_data_subset=variable_data[lat_inds,lon_inds]

    #extraigo mes (climatologia mensual) o meses (climatologia trimestral)
    if (climatologia_tipo=="mensual"):
        meses=str(data_list_climatologia[indice_list].name)[-2]+str(data_list_climatologia[indice_list].name)[-1]
        if (meses=="01"):
            mes="Enero"
        if (meses=="02"):
            mes="Febrero"
        if (meses=="03"):
            mes="Marzo"
        if (meses=="04"):
            mes="Abril"
        if (meses=="05"):
            mes="Mayo"
        if (meses=="06"):
            mes="Junio"
        if (meses=="07"):
            mes="Julio"
        if (meses=="08"):
            mes="Agosto"
        if (meses=="09"):
            mes="Septiembre"
        if (meses=="10"):
            mes="Octubre"
        if (meses=="11"):
            mes="Noviembre"
        if (meses=="12"):
            mes="Diciembre"
            
    if (climatologia_tipo=="trimestral"):
        meses=str(data_list_climatologia[indice_list].name)[-8:-1]+str(data_list_climatologia[indice_list].name)[-1]
        if (meses=="03-04-05"):
            mes="MAM"
        if (meses=="06-07-08"):
            mes="JJA"
        if (meses=="09-10-11"):
            mes="SON"
        if (meses=="12-01-02"):
            mes="DEF"
            
    #ploteo
    if (orientacion=="H"):
        fig1 = plt.figure(figsize=[9,5],dpi=200) #horizontal region 1
    if (orientacion=="V"):
        fig1 = plt.figure(figsize=[7.5,7.5],dpi=200) #vertical sudamerica
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

    if (paleta_color=="matter"):
        variable_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(valor_minimo, valor_maximo, delta_valor),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': variable_data_subset.units},
                   cmap=cmocean.cm.matter)
        
    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.5',linewidth=0.7,alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none', 
                  edgecolor='0.4',alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xticklabels(np.arange(xticks_min,xticks_max)[::espacio_entre_lat_lon])
    plt.xticks(np.arange(xticks_min,xticks_max)[::espacio_entre_lat_lon])
    ax.set_xlabel("Longitud")

    ax.set_yticklabels(np.arange(yticks_min,yticks_max)[::espacio_entre_lat_lon])
    plt.yticks(np.arange(yticks_min,yticks_max)[::espacio_entre_lat_lon])
    ax.set_ylabel("Latitud")

    if (grid==True):
        plt.grid(linestyle="--", alpha=0.3)

    plt.title(variable+" "+region+" "+mes)
    #plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes)
    plt.show()



#%%
"""
Grafico climatologias
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

#ploteo para sudamerica 
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_list,i,"cldamt media mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",8,"V")

import matplotlib.pyplot as plt
plt.close(fig="all")
    
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_list,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_list,i,"cldamt media trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_list,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",8,"V")

plt.close(fig="all")

#ploteo para region 1
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_list,i,"cldamt media mensual (1984-2016)","mensual",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",4,"H")

plt.close(fig="all")

for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_list,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",4,"H")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_list,i,"cldamt media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",4,"H")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_list,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",4,"H")

#ploteo para region 2
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_list,i,"cldamt media mensual (1984-2016)","mensual",-32,-22,-64,-53,"%",0,101,5,-63,-53,-31,-22,True,"Región 2","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",2,"H")

plt.close(fig="all")

for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_list,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-32,-22,-64,-53,"%",0,26,2,-63,-53,-31,-22,True,"Región 2","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",2,"H")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_list,i,"cldamt media trimestral (1984-2016)","trimestral",-32,-22,-64,-53,"%",0,101,5,-63,-53,-31,-22,True,"Región 2","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",2,"H")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_list,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-32,-22,-64,-53,"%",0,26,2,-63,-53,-31,-22,True,"Región 2","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",2,"H")

#%%
"""
veo como seleccionar region con shape
https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan 
https://gis.stackexchange.com/questions/289775/masking-netcdf-data-using-shapefile-xarray-geopandas
"""
import xarray as xr
import rioxarray as rio
from shapely.geometry import mapping, Polygon
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
US=world[world['name']=='United States of America']  # conus, AK, HI, PR
US.geometry.apply(mapping)
veo=provincias[19].exterior.coords[:]
veo=Polygon(veo)
veoveo=MultiPolygon([veo])

veo=IGN.geometry[19]

veoveo=IGN.geometry[19]
IGN.geometry[0].apply(mapping)

#geometries=provincias[19].exterior.coords[:] #esto extrae coordenadas del borde de la provincia. Lo hace como lista, lo puedo pasar a array. 
#geometries=Polygon(geometries)

#geometries=MultiPolygon(geometries)
#xds = media_mensual_cldamt_list[0].rio.write_crs("epsg:4326", inplace=True)
xds = media_mensual_cldamt_list[0]
xds=xds.sel(lat=slice(-32.5,-22.5),lon=slice(360-64.5,360-53.5))
xds=xds.rio.write_crs("epsg:4326",inplace=True)


xds=media_mensual_cldamt_list[0]
xds_veo=xds.rio.write_crs(IGN.crs)

clipped = xds_veo.rio.clip(IGN.geometry.apply(mapping)[19],IGN.crs,drop=False, invert=False)
clipped = xds_veo.rio.clip(IGN.geometry[19][0],IGN.crs)#,drop=False, invert=False)


IGN.geometry.apply(mapping)[19]
type(IGN.geometry.apply(mapping)[19])
veo=IGN.geometry[19][0]

float(veo.wkt)

xds = media_mensual_cldamt_list[0].rio.write_crs("epsg:4326",inplace=True)

newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[IGN.geometry[19][0]])  #esto esta bien, queda igual que tutorial la forma del poligono que recorta https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan

geometria=newconus.geometry.apply(mapping) #lon lat lon lat lon lat

newconus.crs
xds=data_list[0].cldamt
xds=xds.rio.write_crs("epsg:4326", inplace=True)
#xds=xds.transpose(transpose_coords=True)
clipped = xds.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=True) #ESTO TENDRIA QUE FUNCIONAR tiene que haber algo mal en el xds
#esto tira el error este: DimensionError: x dimension not found. 'set_spatial_dims()' can address this.
#voy a cambiar los nombres de lon y lat por x y o algo asi-

#The issue you are facing is that rioxarray expects your spatial dimensions and coordinate to have the same name. I would recommend using the rename methods of in xarray to rename the dimensions and coordinates so they are both longitide and latitude or x and y.

xds=xds.rio.set_spatial_dims("x", "y", inplace=True) #no funciona pruebo otra cosa

###############################ESTAS LINEAS PARECEN FUNCIONAR##########################
newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[IGN.geometry[19][0]])  #esto esta bien, queda igual que tutorial la forma del poligono que recorta https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan
xds=data_list[0].cldamt
xds=xds.swap_dims({"lon": "x"})
xds=xds.swap_dims({"lat": "y"})
xds=xds.rio.write_crs("epsg:4326", inplace=True)
clipped = xds.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False)

#The problem here is that polygon in the geodataframe uses longitudes -180 to 180 whereas my netcdf data uses longitudes 0 to 360. Modifico el xarray a ver...

#xds=xds.assign_coords(lon=(xds.lon - 180)) #fijarme si es asi o es que le tengo que sumar al polygon, creo que es la segunda opcion!! porque estoy siemrpe trabajando de 0 a 360 no? bueno ver desde aca que tal 
#xds.lon

#me genera otras coordenadas x y en el clipped respecto al xds. Busco cambiar los nombres POSTA de las coordenadas. 
#no es esto, lo que esta raro es la latitud!!! ver en grafico como queda. Quizas estan invertidas lon y lats? 
coords=np.array(IGN.geometry[19][0].exterior.coords)
coords[:,0]=coords[:,0]+360 #ver si sumo 180 o 360
coords[:,1]=coords[:,1]+90
newpoly=Polygon(coords)
newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[newpoly])

coords=np.array(IGN.geometry[19][0].exterior.coords)
#invierto coordenadas a ver si se soluciona
veo=coords.copy()
veo[:,1]=coords[:,0]
veo[:,0]=coords[:1]
newpoly=Polygon(veo)
newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[newpoly])
#no funciona

#seguir guiandome de aca: https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan
######################################################################################
    #cargo lats y lons para armar xarray
 #   lats=data_list[0][variable].mean("time", keep_attrs=True)["lat"][:].values
 #   lons=data_list[0][variable].mean("time", keep_attrs=True)["lon"][:].values
 #   coords=[("lat", lats),("lon", lons)]

lats=xds[0]["lat"][:].values
lons=xds[0]["lon"][:].values
coords=[("y",lats),("x",lons)]

xds=xr.DataArray(data_list[0].cldamt.values[0], coords=coords)

#####GRAFICO
plottime='1983-07-15'
fig=plt.figure(figsize=(10,10))
ax1=fig.add_subplot(211,projection=ccrs.PlateCarree())
clipped.sel(time=plottime).plot(ax=ax1)
newconus.boundary.plot(ax=ax1,color='black')
plt.title("clipped pr invert=True")


fig=plt.figure(figsize=(10,10))
ax1=fig.add_subplot(211,projection=ccrs.PlateCarree())
clipped.plot(ax=ax1)
newconus.boundary.plot(ax=ax1,color='black')
plt.title("clipped pr invert=True")


##########################################
clipped = xds.rio.clip(newconus.geometry.apply(mapping)[0],newconus.crs,drop=False, invert=False) #ESTO TENDRIA QUE FUNCIONAR tiene que haber algo mal en el xds


clipped=xds.rio.clip(int(float(str(geometria[0]["coordinates"][0]))),newconus.crs,drop=False, invert=False) #ESTO TENDRIA QUE FUNCIONAR tiene que haber algo mal en el xds
plt.plot(clipped)
clipped = xds.rio.clip(int(float(str(geometria[0]["coordinates"][0]))),newconus.crs,drop=False, invert=False) #ESTO TENDRIA QUE FUNCIONAR tiene que haber algo mal en el xds

clipped = xds_veo.rio.clip(newconus.geometry,newconus.crs,drop=False, invert=False)

xds=xds.rio.set_spatial_dims('lon', 'lat', inplace=True)
obs_dataset_full.rio.set_spatial_dims('x', 'y', inplace=True)


row=next(IGN.iterrows())[1]
row=IGN[IGN.gid==19]
coords=np.array(row['geometry'][0].exterior.coords)
coords[:,0]=coords[:,0]+360.
newpoly=Polygon(coords)
newconus = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[newpoly])

clipped = xds_veo.rio.clip(IGN.geometry.apply(mapping)[19],IGN.crs,drop=False, invert=False)

#ejemplo
world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
conus=world[world['name']=='United States of America']  # conus, AK, HI, PR
conus.geometry

modDS=xr.open_dataset("/home/nadia/Descargas/regrid_pr_Amon_ACCESS1-0_historical_r1i1p1_190601-200512.nc")
pr=modDS.pr
pr=pr.rio.write_crs("epsg:4326", inplace=True)

clipped = pr.rio.clip(IGN.geometry.apply(mapping)[19],IGN.crs,drop=False, invert=False)



#### ESTO ES LO QUE FUNCIONA #####

newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[IGN.geometry[19][0]])  #esto esta bien, queda igual que tutorial la forma del poligono que recorta https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan
xds=data_list[0].cldamt
xds=xds.swap_dims({"lon": "x"})
xds=xds.swap_dims({"lat": "y"})

#clipped = xds.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False)

coords=np.array(IGN.geometry[19][0].exterior.coords)
coords[:,0]=coords[:,0]+360 #ver si sumo 180 o 360
newpoly=Polygon(coords)
newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[newpoly])

lats=xds[0]["lat"][:].values
lons=xds[0]["lon"][:].values
coords=[("y",lats),("x",lons)]

xds=xr.DataArray(data_list[0].cldamt.values[0], coords=coords)
xds=xds.rio.write_crs("epsg:4326", inplace=True)
clipped = xds.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False)

#Pasa que con el tamaño de las grillas quedan muchos lugares en blanco en la provincia. 
#Veo que incluya un area un poco mas grande que corrientes 

## Voy a buscar re grid el xarray para que en lugar de ser de 1ºx1º sea de 0.1ºx0.1º asi queda mejor adentro de corrientes cuando lo clipeo. Lo hago con el metodo de interpolacion del vecino mas cercano https://stackoverflow.com/questions/49973049/how-to-re-gridding-xarray-from-higher-to-lower-resolution-using-idw 

ynuevo=np.linspace(-89.5, 89.5, 3600)
xnuevo=np.linspace(0.5, 359.5, 7200)
xds_2=xds.reindex(y=ynuevo,x=xnuevo, method="nearest") 
xds_3=xds.interp(y=ynuevo,x=xnuevo,method="linear")#ESTA ES LA QUE VA http://xarray.pydata.org/en/stable/generated/xarray.DataArray.interp.html 


clipped = xds_2.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False)
clipped = xds_3.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False) #ESTA ES LA QUE VA

#selecciono region
lats=clipped["y"][:]
lons=clipped["x"][:]
lat_lims=[-31,-27]
lon_lims=[300,305] #lean 360-64 (64 O) 360-31 (31 O) 
lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
clipped_subset=clipped[lat_inds,lon_inds]



fig=plt.figure(figsize=(10,10))
ax1=fig.add_subplot(projection=ccrs.PlateCarree(central_longitude=0))
ax1.coastlines(color='0.3')

clipped_subset.plot.contourf(ax=ax1,
                   levels=np.arange(0, 100, 5),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   #cbar_kwargs={'label': variable_data_subset.units},
                   #cmap=cmocean.cm.rain)
                   )
newconus.boundary.plot(ax=ax1,color='black',transform=ccrs.PlateCarree())

ax1.set_xticklabels(np.arange(-60,-55)[::1])
plt.xticks(np.arange(-60,-55)[::1])
ax1.set_xlabel("Longitud")

ax1.set_yticklabels(np.arange(-31,-26)[::1])
plt.yticks(np.arange(-31,-26)[::1])
ax1.set_ylabel("Latitud")


plt.title("clipped pr invert=True")



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
