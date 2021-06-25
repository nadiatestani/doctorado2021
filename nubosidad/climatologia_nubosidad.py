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
Antes de hacer nada, lo primero que hago es generar en cada archivo xarray tres nuevas variables que van a representar el porcentaje de cldamt de nubes bajas medias y altas, siendo estos las sumas de
cldamt_types : BAJAS:[0:6] MEDIAS:[6:12] ALTAS:[12:19]
"""
import numpy as np

for i in range(0,len(data_list)):
    bajas=data_list[i]["cldamt_types"].values[0][0]+data_list[i]["cldamt_types"].values[0][1]+data_list[i]["cldamt_types"].values[0][2]+data_list[i]["cldamt_types"].values[0][3]+data_list[i]["cldamt_types"].values[0][4]+data_list[i]["cldamt_types"].values[0][5]
    bajas_3D=np.reshape(bajas,(1,180,360))
    bajas_array=xr.DataArray(data=bajas_3D,dims=["time","lat","lon"])
    bajas_array=bajas_array.assign_coords({"lon":data_list[i].lon.values ,"lat":data_list[i].lat.values ,"time":data_list[i].time.values })
    bajas_array.attrs["units"]= "percent"
    
    medias=data_list[i]["cldamt_types"].values[0][6]+data_list[i]["cldamt_types"].values[0][7]+data_list[i]["cldamt_types"].values[0][8]+data_list[i]["cldamt_types"].values[0][9]+data_list[i]["cldamt_types"].values[0][10]+data_list[i]["cldamt_types"].values[0][11]
    medias_3D=np.reshape(medias,(1,180,360))
    medias_array=xr.DataArray(data=medias_3D,dims=["time","lat","lon"])
    medias_array=medias_array.assign_coords({"lon":data_list[i].lon.values ,"lat":data_list[i].lat.values ,"time":data_list[i].time.values })
    medias_array.attrs["units"]= "percent"
    
    altas=data_list[i]["cldamt_types"].values[0][12]+data_list[i]["cldamt_types"].values[0][13]+data_list[i]["cldamt_types"].values[0][14]+data_list[i]["cldamt_types"].values[0][15]+data_list[i]["cldamt_types"].values[0][16]+data_list[i]["cldamt_types"].values[0][17]
    altas_3D=np.reshape(altas,(1,180,360))
    altas_array=xr.DataArray(data=altas_3D,dims=["time","lat","lon"])
    altas_array=altas_array.assign_coords({"lon":data_list[i].lon.values ,"lat":data_list[i].lat.values ,"time":data_list[i].time.values })
    altas_array.attrs["units"]= "percent"
    
    data_list[i]=data_list[i].assign(cldamt_bajas=bajas_array)
    data_list[i]=data_list[i].assign(cldamt_medias=medias_array)
    data_list[i]=data_list[i].assign(cldamt_altas=altas_array)



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

#extiendo a cldamt_bajas, cldamt_medias, cldamt_altas

#armo lista con medias mensuales de cldamt_bajas
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
media_mensual_cldamt_bajas_list=[None]*12
for i in range(0,12):
    media_mensual_cldamt_bajas_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_bajas",meses[i],"%")[0]

#armo lista con desvios mensuales de cldamt_bajas
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
desvio_mensual_cldamt_bajas_list=[None]*12
for i in range(0,12):
    desvio_mensual_cldamt_bajas_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_bajas",meses[i],"%")[1]

#armo lista con medias mensuales de cldamt_medias
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
media_mensual_cldamt_medias_list=[None]*12
for i in range(0,12):
    media_mensual_cldamt_medias_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_medias",meses[i],"%")[0]

#armo lista con desvios mensuales de cldamt_medias
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
desvio_mensual_cldamt_medias_list=[None]*12
for i in range(0,12):
    desvio_mensual_cldamt_medias_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_medias",meses[i],"%")[1]

#armo lista con medias mensuales de cldamt_altas
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
media_mensual_cldamt_altas_list=[None]*12
for i in range(0,12):
    media_mensual_cldamt_altas_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_altas",meses[i],"%")[0]

#armo lista con desvios mensuales de cldamt_altas
meses=["01","02","03","04","05","06","07","08","09","10","11","12"]
desvio_mensual_cldamt_altas_list=[None]*12
for i in range(0,12):
    desvio_mensual_cldamt_altas_list[i]=media_desvio_mensual(data_list_1984_2016,"cldamt_altas",meses[i],"%")[1]


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

#extiendo a cldamt bajas medias y altas

#armo lista con medias trimestrales cldamt_bajas
media_trimestral_cldamt_bajas_list=[None]*4
media_trimestral_cldamt_bajas_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","03","04","05","%")[0]
media_trimestral_cldamt_bajas_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","06","07","08","%")[0]
media_trimestral_cldamt_bajas_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","09","10","11","%")[0]
media_trimestral_cldamt_bajas_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","12","01","02","%")[0]

#armo lista con desvios trimestrales cldamt_bajas
desvio_trimestral_cldamt_bajas_list=[None]*4
desvio_trimestral_cldamt_bajas_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","03","04","05","%")[1]
desvio_trimestral_cldamt_bajas_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","06","07","08","%")[1]
desvio_trimestral_cldamt_bajas_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","09","10","11","%")[1]
desvio_trimestral_cldamt_bajas_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_bajas","12","01","02","%")[1]

#armo lista con medias trimestrales cldamt_medias
media_trimestral_cldamt_medias_list=[None]*4
media_trimestral_cldamt_medias_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","03","04","05","%")[0]
media_trimestral_cldamt_medias_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","06","07","08","%")[0]
media_trimestral_cldamt_medias_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","09","10","11","%")[0]
media_trimestral_cldamt_medias_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","12","01","02","%")[0]

#armo lista con desvios trimestrales cldamt_medias
desvio_trimestral_cldamt_medias_list=[None]*4
desvio_trimestral_cldamt_medias_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","03","04","05","%")[1]
desvio_trimestral_cldamt_medias_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","06","07","08","%")[1]
desvio_trimestral_cldamt_medias_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","09","10","11","%")[1]
desvio_trimestral_cldamt_medias_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_medias","12","01","02","%")[1]

#armo lista con medias trimestrales cldamt_altas
media_trimestral_cldamt_altas_list=[None]*4
media_trimestral_cldamt_altas_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","03","04","05","%")[0]
media_trimestral_cldamt_altas_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","06","07","08","%")[0]
media_trimestral_cldamt_altas_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","09","10","11","%")[0]
media_trimestral_cldamt_altas_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","12","01","02","%")[0]

#armo lista con desvios trimestrales cldamt_altas
desvio_trimestral_cldamt_altas_list=[None]*4
desvio_trimestral_cldamt_altas_list[0]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","03","04","05","%")[1]
desvio_trimestral_cldamt_altas_list[1]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","06","07","08","%")[1]
desvio_trimestral_cldamt_altas_list[2]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","09","10","11","%")[1]
desvio_trimestral_cldamt_altas_list[3]=media_desvio_trimestral(data_list_12_1983_11_2016,"cldamt_altas","12","01","02","%")[1]



#%%
"""
Recorto los xarrays de la climatologia para el shape de Corrientes
Armo funcion que selecciona region de xarray con shape
https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan 
https://gis.stackexchange.com/questions/289775/masking-netcdf-data-using-shapefile-xarray-geopandas
"""

import numpy as np
import xarray as xr
import rioxarray as rio
from shapely.geometry import mapping, Polygon
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def clip_xarray_con_shape(coords,xarray_lista, numero_elemento_lista,climatologia_o_original,variable_original):
    """
    

    Parameters
    ----------
    coords : np.array
        Coordenadas del borde del shape que se quiere recortar. Por ejemplo: IGN=geopandas.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp"); coords=np.array(IGN.geometry[19][0].exterior.coords)
    xarray_lista: list
        lista con arrays que se quieren recortar. Es decir: media_mensual_cldamt_list[numero_elemento_lista] 
    numero_elemento_lista: float
        numero del elemento de la lista
    climatologia_o_original: str
        climatologia: para trabajar con xarrays de climatologia que tienen una unica variable
        original: para trabajar con xarrays originales
    variable_original: str
        variable que se esta trabajando, va a usarse solo en la funcion si se carga la lista con xarrays originales
    Returns
    -------
    List con primer elemento: xarray recortado, segundo elemento: shape que se recorto

    """
    coords[:,0]=coords[:,0]+360 #sumo 360 a las longitudes para que tenga la misma referencia que los xarrays
    newpoly=Polygon(coords) #armo poligono
    newconus=gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[newpoly]) #georeferencio el poligono

    #cargo xarray para recortar cambio el nombre de las dimensiones a x,y 
    if (climatologia_o_original=="climatologia"):
        xds=xarray_lista[numero_elemento_lista]
        xds_nombre=xds.name
        
    if (climatologia_o_original=="original"):
        xds=xarray_lista[numero_elemento_lista][variable_original]
        xds_nombre=xds.name+" "+str(xds["time"].values[0])[0:10]
    
    xds=xds.swap_dims({"lon": "x"})
    xds=xds.swap_dims({"lat": "y"})

    lats=xds["lat"][:].values
    lons=xds["lon"][:].values
    coords_2=[("y",lats),("x",lons)]
    
    if (climatologia_o_original=="climatologia"):
        xds=xr.DataArray(xarray_lista[numero_elemento_lista].values, coords=coords_2)
    if (climatologia_o_original=="original"):
        xds=xr.DataArray(xarray_lista[numero_elemento_lista][variable_original].values[0], coords=coords_2)
    
    xds=xds.rio.write_crs("epsg:4326", inplace=True) #georeferencio
    
    #lo regrillo para que pase de ser 1ºx1º a 0.05ºx0.05º asi se recorta bien. Lo hago con el metodo de interpolacion lineal. el regrillado lo hago solo en un entorno a corrientes para ahorrar memoria
    #ynuevo=np.linspace(-89.5, 89.5, 3600)
    #xnuevo=np.linspace(0.5, 359.5, 7200)
    ynuevo=np.linspace(-31.5, -25.5, 60)
    xnuevo=np.linspace(280.5, 315.5, 350)
    xds=xds.interp(y=ynuevo,x=xnuevo,method="linear")# http://xarray.pydata.org/en/stable/generated/xarray.DataArray.interp.html 
    
    clipped = xds.rio.clip(newconus.geometry.apply(mapping),newconus.crs,drop=False, invert=False) #lo recorto
    clipped.name=xds_nombre
    #clipped.attrs['units']=xds.attrs['units']
    #clipped.attrs['mes']=xds.attrs['mes']
    return(clipped,newconus)

#%%
"""
Corro esta funcion y armo lista con los estadisticos (media y desvio) recortados en la provincia de Corrientes mensual y trimestral y otra lista con todos los dias clippeados de cldamt
"""

IGN=gpd.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")

#armo lista con variable original, medias y desvios mensuales de cldamt

cldamt_list_clipped=[None]*len(data_list)
for i in range(0,len(data_list)):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    cldamt_list_clipped[i]=clip_xarray_con_shape(coords, data_list, i,"original","cldamt")[0]


media_mensual_cldamt_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_mensual_cldamt_list_clipped[i]=clip_xarray_con_shape(coords, media_mensual_cldamt_list, i,"climatologia","cldamt")[0]

desvio_mensual_cldamt_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_mensual_cldamt_list_clipped[i]=clip_xarray_con_shape(coords, desvio_mensual_cldamt_list, i,"climatologia","cldamt")[0]

#armo lista con medias y desvios trimestrales de cldamt

media_trimestral_cldamt_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_trimestral_cldamt_list_clipped[i]=clip_xarray_con_shape(coords, media_trimestral_cldamt_list, i,"climatologia","cldamt")[0]

desvio_trimestral_cldamt_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_trimestral_cldamt_list_clipped[i]=clip_xarray_con_shape(coords, desvio_trimestral_cldamt_list, i,"climatologia","cldamt")[0]


#extiendo a cldamt_bajas, cldamt_medias, cldamt_altas

###########bajas
cldamt_bajas_list_clipped=[None]*len(data_list)
for i in range(0,len(data_list)):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    cldamt_bajas_list_clipped[i]=clip_xarray_con_shape(coords, data_list, i,"original","cldamt_bajas")[0]


media_mensual_cldamt_bajas_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_mensual_cldamt_bajas_list_clipped[i]=clip_xarray_con_shape(coords, media_mensual_cldamt_bajas_list, i,"climatologia","cldamt_bajas")[0]

desvio_mensual_cldamt_bajas_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_mensual_cldamt_bajas_list_clipped[i]=clip_xarray_con_shape(coords, desvio_mensual_cldamt_bajas_list, i,"climatologia","cldamt_bajas")[0]

#armo lista con medias y desvios trimestrales de cldamt_bajas

media_trimestral_cldamt_bajas_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_trimestral_cldamt_bajas_list_clipped[i]=clip_xarray_con_shape(coords, media_trimestral_cldamt_bajas_list, i,"climatologia","cldamt_bajas")[0]

desvio_trimestral_cldamt_bajas_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_trimestral_cldamt_bajas_list_clipped[i]=clip_xarray_con_shape(coords, desvio_trimestral_cldamt_bajas_list, i,"climatologia","cldamt_bajas")[0]

###########medias
cldamt_medias_list_clipped=[None]*len(data_list)
for i in range(0,len(data_list)):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    cldamt_medias_list_clipped[i]=clip_xarray_con_shape(coords, data_list, i,"original","cldamt_medias")[0]


media_mensual_cldamt_medias_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_mensual_cldamt_medias_list_clipped[i]=clip_xarray_con_shape(coords, media_mensual_cldamt_medias_list, i,"climatologia","cldamt_medias")[0]

desvio_mensual_cldamt_medias_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_mensual_cldamt_medias_list_clipped[i]=clip_xarray_con_shape(coords, desvio_mensual_cldamt_medias_list, i,"climatologia","cldamt_medias")[0]

#armo lista con medias y desvios trimestrales de cldamt_medias

media_trimestral_cldamt_medias_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_trimestral_cldamt_medias_list_clipped[i]=clip_xarray_con_shape(coords, media_trimestral_cldamt_medias_list, i,"climatologia","cldamt_medias")[0]

desvio_trimestral_cldamt_medias_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_trimestral_cldamt_medias_list_clipped[i]=clip_xarray_con_shape(coords, desvio_trimestral_cldamt_medias_list, i,"climatologia","cldamt_medias")[0]

###########altas
cldamt_altas_list_clipped=[None]*len(data_list)
for i in range(0,len(data_list)):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    cldamt_altas_list_clipped[i]=clip_xarray_con_shape(coords, data_list, i,"original","cldamt_altas")[0]


media_mensual_cldamt_altas_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_mensual_cldamt_altas_list_clipped[i]=clip_xarray_con_shape(coords, media_mensual_cldamt_altas_list, i,"climatologia","cldamt_altas")[0]

desvio_mensual_cldamt_altas_list_clipped=[None]*12
for i in range(0,12):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_mensual_cldamt_altas_list_clipped[i]=clip_xarray_con_shape(coords, desvio_mensual_cldamt_altas_list, i,"climatologia","cldamt_altas")[0]

#armo lista con medias y desvios trimestrales de cldamt_altas

media_trimestral_cldamt_altas_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    media_trimestral_cldamt_altas_list_clipped[i]=clip_xarray_con_shape(coords, media_trimestral_cldamt_altas_list, i,"climatologia","cldamt_altas")[0]

desvio_trimestral_cldamt_altas_list_clipped=[None]*4
for i in range(0,4):
    coords=np.array(IGN.geometry[19][0].exterior.coords)
    desvio_trimestral_cldamt_altas_list_clipped[i]=clip_xarray_con_shape(coords, desvio_trimestral_cldamt_altas_list, i,"climatologia","cldamt_altas")[0]



#%%
"""
Defino funcion que grafica climatologia mensual de alguna variable definida previamente en una determinada region (region cuadrada)
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
        #######ver estas lineas###############
    facecolors = np.ma.array(variable_data_subset, mask=np.isnan(variable_data_subset))
    ax.set_array(facecolors)
        #######ver estas lineas###############
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
Defino funcion que grafica climatologia mensual de alguna variable definida previamente en una determinada region (region con shape)
"""

def grafico_campos_climatologia_nubosidad_clip(paises,provincias,data_list_climatologia,indice_list,variable,climatologia_tipo,lat_min,lat_max,lon_min,lon_max,unidades_nombre,valor_minimo, valor_maximo, delta_valor,xticks_min,xticks_max, yticks_min, yticks_max,grid,region, ruta_salida, paleta_color,espacio_entre_lat_lon,orientacion):
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
    #variable_data.attrs["units"]=unidades_nombre #cambio el nombre de la unidad

    #selecciono region
    lats=variable_data["y"][:]
    lons=variable_data["x"][:]
    lat_lims=[lat_min,lat_max]
    lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
    lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
    lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
    variable_data_subset=variable_data[lat_inds,lon_inds]

    #extraigo mes (climatologia mensual) o meses (climatologia trimestral)
    if (climatologia_tipo=="mensual"):
            meses=str(data_list_climatologia[indice_list].name)[-2]+str(data_list_climatologia[indice_list].name)[-1]
            if (meses=="01"):
                mes1="Enero"
            if (meses=="02"):
                mes1="Febrero"
            if (meses=="03"):
                mes1="Marzo"
            if (meses=="04"):
                mes1="Abril"
            if (meses=="05"):
                mes1="Mayo"
            if (meses=="06"):
                mes1="Junio"
            if (meses=="07"):
                mes1="Julio"
            if (meses=="08"):
                mes1="Agosto"
            if (meses=="09"):
                mes1="Septiembre"
            if (meses=="10"):
                mes1="Octubre"
            if (meses=="11"):
                mes1="Noviembre"
            if (meses=="12"):
                mes1="Diciembre"

    if (climatologia_tipo=="trimestral"):
        meses=str(data_list_climatologia[indice_list].name)[-8:-1]+str(data_list_climatologia[indice_list].name)[-1]
        if (meses=="03-04-05"):
            mes1="MAM"
        if (meses=="06-07-08"):
            mes1="JJA"
        if (meses=="09-10-11"):
            mes1="SON"
        if (meses=="12-01-02"):
            mes1="DEF"
    
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
                   cbar_kwargs={'label': unidades_nombre},
                   cmap=cmocean.cm.rain)
    
    if (paleta_color=="curl"):
        variable_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(valor_minimo, valor_maximo, delta_valor),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': unidades_nombre},
                   cmap=cmocean.cm.curl_r)

    if (paleta_color=="matter"):
        variable_data_subset.plot.contourf(ax=ax,
                   levels=np.arange(valor_minimo, valor_maximo, delta_valor),
                   extend='neither',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': unidades_nombre},
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

    #plt.title(variable+" "+region)
    plt.title(variable+" "+region+" "+mes1)
    ##plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes1)
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
paises=MultiPolygon([df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == 'Brazil']['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == 'Paraguay']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Uruguay']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Bolivia']['geometry'].values[0],
                     df.loc[df['ADMIN'] == 'Chile']['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Colombia"]['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Ecuador"]['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Venezuela"]['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Guyana"]['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Suriname"]['geometry'].values[0],
                     df.loc[df['ADMIN'] == "Panama"]['geometry'].values[0][0],
                     df.loc[df['ADMIN'] == "Costa Rica"]['geometry'].values[0][0]]) #los paso a multipolygon para poder graficarlos

#cargo shape con provincias de argentina con datos del IGN 
#descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN=geopandas.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias=[None]*24
for i in range(0,24):
    provincias[i]=IGN["geometry"][i]
provincias=MultiPolygon(provincias) #paso a multipolygon para poder ponerlo en mapa

#%%
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

#PLOTEO PARA CORRIENTES
#plt.close(fig="all")

#paleta1
#for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_list_clipped,i,"cldamt media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,101,5,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",1,"H")

plt.close(fig="all")

#paleta2
for i in range(0,12):
    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_list_clipped,i,"cldamt media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",40,70,1,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",1,"H")


#plt.close(fig="all")

#for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_list_clipped,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,26,2,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",1,"H")

#paleta2
plt.close(fig="all")

for i in range(0,12):
    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_list_clipped,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,14,1,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",1,"H")

#plt.close(fig="all")

#for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_list,i,"cldamt media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",4,"H")

#paleta2
plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_trimestral_cldamt_list_clipped,i,"cldamt media trimestral (1984-2016)","trimestral",-31,-26,-60,-55,"%",40,70,1,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",1,"H")


#plt.close(fig="all")

#for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_list,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",4,"H")

#paleta2
plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_trimestral_cldamt_list_clipped,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-31,-26,-60,-55,"%",0,14,1,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",1,"H")


#%% 
#extiendo a cldamt bajas medias y altas
#ploteo para sudamerica 
#bajas
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_bajas_list,i,"cldamt nubes bajas media mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","rain",8,"V")
#ver que quedan espacios en blanco. 

import matplotlib.pyplot as plt
plt.close(fig="all")
    
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_bajas_list,i,"cldamt nubes bajas desvío estándar mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","matter",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_bajas_list,i,"cldamt nubes bajas media trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","rain",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_bajas_list,i,"cldamt nubes bajas desvío estándar trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","matter",8,"V")

plt.close(fig="all")

###############################################################################################################################
#medias
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_medias_list,i,"cldamt nubes medias media mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","rain",8,"V")
#ver que quedan espacios en blanco. 

import matplotlib.pyplot as plt
plt.close(fig="all")
    
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_medias_list,i,"cldamt nubes medias desvío estándar mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","matter",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_medias_list,i,"cldamt nubes medias media trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","rain",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_medias_list,i,"cldamt nubes medias desvío estándar trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","matter",8,"V")

plt.close(fig="all")

###############################################################################################################################
#altas
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_mensual_cldamt_altas_list,i,"cldamt nubes altas media mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","rain",8,"V")
#ver que quedan espacios en blanco. 

import matplotlib.pyplot as plt
plt.close(fig="all")
    
for i in range(0,12):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_mensual_cldamt_altas_list,i,"cldamt nubes altas desvío estándar mensual (1984-2016)","mensual",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","matter",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_altas_list,i,"cldamt nubes altas media trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","rain",8,"V")

plt.close(fig="all")

for i in range(0,4):
    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_altas_list,i,"cldamt nubes altas desvío estándar trimestral (1984-2016)","trimestral",-60,15,-90,-30,"%",0,26,2,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","matter",8,"V")

plt.close(fig="all")



#%%

"""
calculo la media mensual de una determinada variable en determinada region 
grafico una serie temporal de todo el periodo, y por mes
"""
#defino funcion que calcula media mensual de una determinada variable en determinada region

def media_espacial(data_list,indice_list,variable,lat_min,lat_max,lon_min,lon_max,clipped):
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
    clipped: True or False
        True: la entrada es la lista de la variable clipeada
        False: la entrada es la lista sin clipear
    Returns float
    -------
    media espacial de la variable seleccionada en la region seleccionada
    """
    import numpy as np
    
    data=data_list[indice_list]
    if (clipped==False):
        variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
        #selecciono region
        lats=variable_data["lat"][:]
        lons=variable_data["lon"][:]
        lat_lims=[lat_min,lat_max]
        lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
        lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
        lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
        variable_data_subset=variable_data[lat_inds,lon_inds]
        
    if (clipped==True):
        variable_data=data
        #selecciono region
        lats=variable_data["y"][:]
        lons=variable_data["x"][:]
        lat_lims=[lat_min,lat_max]
        lon_lims=[360+lon_min,360+lon_max] #lean 360-64 (64 O) 360-31 (31 O) 
        lat_inds=np.where((lats>lat_lims[0]) & (lats<lat_lims[1]))[0]
        lon_inds=np.where((lons>lon_lims[0]) & (lons<lon_lims[1]))[0]
        variable_data_subset=variable_data[lat_inds,lon_inds]
    
    #calculo media espacial del subset de la variable data
    media_espacial=np.nanmean(variable_data_subset.values,axis=(0,1))
    return(media_espacial)

#armo funcion que devuelve data frame donde primera columna sea la fecha y la segunda columna sea la media espacial
def media_espacial_df(data_list,variable,lat_min,lat_max,lon_min,lon_max,clipped):
    """
    

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, cargar la lista para modificar
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
    clipped: True or False
        True: la entrada es la lista de la variable clipeada
        False: la entrada es la lista sin clipear

    Returns
    -------
    Data frame con primer columna fecha segunda columna la media espacial de la variable elegida para esa fecha

    """
    media_espacial_df=pd.DataFrame(columns=["fecha","Media_espacial_"+variable])
    if (clipped==False):
        for i in range(0,len(data_list)):
            #extraigo mes
            mes=str(data_list[i]["time"].values[0])[5:7]
            #extraigo anio
            anio=str(data_list[i]["time"].values[0])[0:4]
            #¢alculo media
            media_espacial_i=media_espacial(data_list,i,variable,lat_min,lat_max,lon_min,lon_max,False)
            media_espacial_df=media_espacial_df.append({"fecha": mes+"-"+anio, "Media_espacial_"+variable: media_espacial_i},ignore_index=True)
    
    if (clipped==True):
        for i in range(0,len(data_list)):
            #extraigo mes
            mes=str(data_list[i].name)[-5:-3]
            #extraigo anio
            anio=str(data_list[i].name)[-10:-6]
            #¢alculo media
            media_espacial_i=media_espacial(data_list,i,variable,lat_min,lat_max,lon_min,lon_max,True)
            media_espacial_df=media_espacial_df.append({"fecha": mes+"-"+anio, "Media_espacial_"+variable: media_espacial_i},ignore_index=True)
            
    media_espacial_df["fecha"]=pd.to_datetime(media_espacial_df["fecha"])
    return(media_espacial_df)

cldamt_media_espacial_df_sudamerica=media_espacial_df(data_list,"cldamt",-60,15,-90,-30,False)
cldamt_media_espacial_df_region1=media_espacial_df(data_list,"cldamt",-39,-16,-64,-31,False)
cldamt_media_espacial_df_region2=media_espacial_df(data_list,"cldamt",-32,-22,-64,-53,False)
cldamt_media_espacial_df_corrientes=media_espacial_df(cldamt_list_clipped,"cldamt",-32,-22,-64,-53,True)


#separo por mes
def media_espacial_mes(media_espacial_df,numero_mes):
        media_espacial_df_mes=media_espacial_df[(pd.DatetimeIndex(media_espacial_df["fecha"]).month)==numero_mes]
        return(media_espacial_df_mes)

#armo una lista con los df de media espacial mensual dde cldamt para cada region donde cada elemento tenga la serie de cada mes

#sudamerica
data_list_cldamt_media_espacial_mensuales_sudamerica=[None]*12
for i in range(0,12):
    data_list_cldamt_media_espacial_mensuales_sudamerica[i]=media_espacial_mes(cldamt_media_espacial_df_sudamerica,i+1)

#region1
data_list_cldamt_media_espacial_mensuales_region1=[None]*12
for i in range(0,12):
    data_list_cldamt_media_espacial_mensuales_region1[i]=media_espacial_mes(cldamt_media_espacial_df_region1,i+1)

#region2
data_list_cldamt_media_espacial_mensuales_region2=[None]*12
for i in range(0,12):
    data_list_cldamt_media_espacial_mensuales_region2[i]=media_espacial_mes(cldamt_media_espacial_df_region2,i+1)

#corrientes
data_list_cldamt_media_espacial_mensuales_corrientes=[None]*12
for i in range(0,12):
    data_list_cldamt_media_espacial_mensuales_corrientes[i]=media_espacial_mes(cldamt_media_espacial_df_corrientes,i+1)


#separo por estacion
def media_espacial_estacion(media_espacial_df,numero_mes1,numero_mes2,numero_mes3):
        media_espacial_df_meses_estacion=media_espacial_df[(((pd.DatetimeIndex(media_espacial_df["fecha"]).month)==numero_mes1) | ((pd.DatetimeIndex(media_espacial_df["fecha"]).month)==numero_mes2) | ((pd.DatetimeIndex(media_espacial_df["fecha"]).month)==numero_mes3))]
        media_espacial_df_meses_estacion.index=range(0,len(media_espacial_df_meses_estacion))
        media_espacial_df_media_estacion=pd.DataFrame(columns=["fecha","Media_espacial_media_estacion"])
        for i in range(0,len(media_espacial_df_meses_estacion),3):
            media_espacial_df_media_estacion=media_espacial_df_media_estacion.append({"fecha": pd.DatetimeIndex(media_espacial_df_meses_estacion["fecha"]).year[i],"Media_espacial_media_estacion": ((media_espacial_df_meses_estacion["Media_espacial_cldamt"][i]+media_espacial_df_meses_estacion["Media_espacial_cldamt"][i+1]+media_espacial_df_meses_estacion["Media_espacial_cldamt"][i+2])/3)}, ignore_index=True)
        return(media_espacial_df_media_estacion)

#armo una lista con los df de media espacial promedio estacional para cada anio dde cldamt para cada region donde cada elemento tenga la serie de cada estacion

#sudamerica
data_list_cldamt_media_espacial_estacional_sudamerica=[None]*4
data_list_cldamt_media_espacial_estacional_sudamerica[0]=media_espacial_estacion(cldamt_media_espacial_df_sudamerica,12,1,2)
data_list_cldamt_media_espacial_estacional_sudamerica[1]=media_espacial_estacion(cldamt_media_espacial_df_sudamerica,3,4,5)
data_list_cldamt_media_espacial_estacional_sudamerica[2]=media_espacial_estacion(cldamt_media_espacial_df_sudamerica,6,7,8)
data_list_cldamt_media_espacial_estacional_sudamerica[3]=media_espacial_estacion(cldamt_media_espacial_df_sudamerica,9,10,11)

#region1
data_list_cldamt_media_espacial_estacional_region1=[None]*4
data_list_cldamt_media_espacial_estacional_region1[0]=media_espacial_estacion(cldamt_media_espacial_df_region1,12,1,2)
data_list_cldamt_media_espacial_estacional_region1[1]=media_espacial_estacion(cldamt_media_espacial_df_region1,3,4,5)
data_list_cldamt_media_espacial_estacional_region1[2]=media_espacial_estacion(cldamt_media_espacial_df_region1,6,7,8)
data_list_cldamt_media_espacial_estacional_region1[3]=media_espacial_estacion(cldamt_media_espacial_df_region1,9,10,11)

#region2
data_list_cldamt_media_espacial_estacional_region2=[None]*4
data_list_cldamt_media_espacial_estacional_region2[0]=media_espacial_estacion(cldamt_media_espacial_df_region2,12,1,2)
data_list_cldamt_media_espacial_estacional_region2[1]=media_espacial_estacion(cldamt_media_espacial_df_region2,3,4,5)
data_list_cldamt_media_espacial_estacional_region2[2]=media_espacial_estacion(cldamt_media_espacial_df_region2,6,7,8)
data_list_cldamt_media_espacial_estacional_region2[3]=media_espacial_estacion(cldamt_media_espacial_df_region2,9,10,11)

#corrientes
data_list_cldamt_media_espacial_estacional_corrientes=[None]*4
data_list_cldamt_media_espacial_estacional_corrientes[0]=media_espacial_estacion(cldamt_media_espacial_df_corrientes,12,1,2)
data_list_cldamt_media_espacial_estacional_corrientes[1]=media_espacial_estacion(cldamt_media_espacial_df_corrientes,3,4,5)
data_list_cldamt_media_espacial_estacional_corrientes[2]=media_espacial_estacion(cldamt_media_espacial_df_corrientes,6,7,8)
data_list_cldamt_media_espacial_estacional_corrientes[3]=media_espacial_estacion(cldamt_media_espacial_df_corrientes,9,10,11)


#%%
"""
grafico las series temporales
"""
#armo funcion para graficar serie completa
#armo funcion para graficar series de cada estacion (en un mismo plot)
#armo funcion para graficar las series por cada anio (en un mismo plot)

#anual
import matplotlib.pyplot as plt
def serie_periodo_completo(data_frame_entrada,variable,region):
    fig1, ax = plt.subplots(figsize=[10,5],dpi=200)
    plt.plot(data_frame_entrada["fecha"],data_frame_entrada["Media_espacial_cldamt"],color="indigo")
    ax.tick_params(axis='x',direction='out',bottom=True,labelrotation=45, labelsize=10,pad=1.5)
    #ax.set_ylim(20,80)
    ax.set_xlabel("Fecha", size=10)
    ax.set_ylabel("cldamt %", size=10)
    ax.grid()
    plt.title(variable+" Media mensual media "+region+ " (serie completa)")
    nombre=variable+"_"+"media_espacial_mensual_"+region+"_"+"(serie completa)"
    plt.savefig("/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"+nombre, dpi=140)
    plt.show

#corro la funcion para graficar la serie de cada region
serie_periodo_completo(cldamt_media_espacial_df_sudamerica,"cldamt","Sudamérica")
serie_periodo_completo(cldamt_media_espacial_df_region1,"cldamt","Región 1")
serie_periodo_completo(cldamt_media_espacial_df_region2,"cldamt","Región 2")
serie_periodo_completo(cldamt_media_espacial_df_corrientes,"cldamt","Corrientes")


#mensual
import matplotlib.pyplot as plt

def serie_mensual(lista,variable,region,ymin,ymax):

    fig1, ax = plt.subplots(4,3,figsize=[12,10],dpi=200) #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1.suptitle(variable+" Media mensual media "+region+ " (meses)",size=18)
    
    meses=[["Enero","Febrero","Marzo"],["Abril","Mayo","Junio"],["Julio","Agosto","Septiembre"],["Octubre","Noviembre","Diciembre"]]
    for j in range(0,4):
        for i in range(0,3):
            ax[j,i].plot(lista[i+3*j]["fecha"],lista[i+3*j]["Media_espacial_cldamt"],color="indigo")
            ax[j,i].tick_params(axis='x',direction='out',bottom=True,labelrotation=45, labelsize=10,pad=1.5)
            ax[j,i].set_ylim(ymin,ymax)
            ax[j,i].set_xlabel("Fecha", size=10)
            ax[j,i].set_ylabel("cldamt %", size=10)
            ax[j,i].grid()
            ax[j,i].set_title(meses[j][i])

    fig1.tight_layout()
    nombre=variable+"_"+"media_espacial_mensual_"+region+"_"+"(meses)"
    plt.savefig("/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"+nombre, dpi=140)
    plt.show

serie_mensual(data_list_cldamt_media_espacial_mensuales_sudamerica,"cldamt","Sudamérica",60,80)
serie_mensual(data_list_cldamt_media_espacial_mensuales_region1,"cldamt","Región 1",50,80)
serie_mensual(data_list_cldamt_media_espacial_mensuales_region2,"cldamt","Región 2",30,80)
serie_mensual(data_list_cldamt_media_espacial_mensuales_corrientes,"cldamt","Corrientes",20,80)

#trimestral

import matplotlib.pyplot as plt

def serie_trimestral(lista,variable,region,ymin,ymax):

    fig1, ax = plt.subplots(2,2,figsize=[12,10],dpi=200) #https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1.suptitle(variable+" Media trimestral media "+region+ " (meses)",size=18)
    
    estacion=[["DEF","MAM"],["JJA","SON"]]
    for j in range(0,2):
        for i in range(0,2):
            ax[j,i].plot(lista[i+2*j]["fecha"],lista[i+2*j]["Media_espacial_media_estacion"],color="indigo")
            ax[j,i].tick_params(axis='x',direction='out',bottom=True,labelrotation=45, labelsize=10,pad=1.5)
            ax[j,i].set_ylim(ymin,ymax)
            ax[j,i].set_xlabel("Fecha", size=10)
            ax[j,i].set_ylabel("cldamt %", size=10)
            ax[j,i].grid()
            ax[j,i].set_title(estacion[j][i])

    fig1.tight_layout()
    nombre=variable+"_"+"media_espacial_trimestral_"+region+"_"+"(estaciones)"
    plt.savefig("/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"+nombre, dpi=140)
    plt.show

serie_trimestral(data_list_cldamt_media_espacial_estacional_sudamerica,"cldamt","Sudamérica",60,80)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region1,"cldamt","Región 1",50,80)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region2,"cldamt","Región 2",40,80)
serie_trimestral(data_list_cldamt_media_espacial_estacional_corrientes,"cldamt","Corrientes",40,80)

#%%
"""
LISTO CON LAS SERIES, extender a tipos de nube todo el analisis. 
"""
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
