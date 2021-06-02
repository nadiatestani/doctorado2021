#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Climatologia de nubosidad con datos de ISCCP-H HGM 

Los datos vienen dados como Promedios mensuales con una resolucion espacial de 1º, estan en formato netCDF CF complaint y se descargaron directamente o con THREDDS. Los detalles de los datos estan en este link: https://www.ncdc.noaa.gov/isccp/isccp-data-access/isccp-basic-data 

Los promedios mensuales HGM surgen de hacer el promedio mensual de los datos ISCCP-H HGH que son los datos del promedio cada 3 horas para cada hora sinoptica. 

Los datos empiezan en el mes 07 del año 1983 y se extienden hasta el mes 06 del año 2017

Para descargarlos use la opcion de descarga con THREEDS 

Para descargar la informacion use la opcion de descarga directa con este comando en la terminal (fuente: https://stackoverflow.com/questions/6827459/can-i-use-wget-to-download-multiple-files-from-linux-terminal)

>wget -r -l1 -A.nc https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access/isccp-basic/hgm/

"""

#%% 
"""
Visualizo la informacion:
    Para eso abro uno de los archivos .nc de nubosidad mensual
    Veo las variables
    Hago un primer grafico exploratorio 
    Uso este tutorial https://carpentrieslab.github.io/python-aos-lesson/02-visualisation/index.html
"""

import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np


nc_ruta="Documentos/Doctorado/datos/nubosidad/ISCCP-H_HGM"
nc_name="ISCCP-Basic.HGM.v01r00.GLOBAL.1983.07.99.9999.GPC.10KM.CS00.EA1.00.nc"
dset=xr.open_dataset(nc_ruta+"/"+nc_name)
print(dset)
print(dset["cldamt"]) #Mean cloud amount
print(dset["cldamt_dist"])
print(dset["pc"])
print(dset["tc"])
print(dset["sigma_tc_time"])
print(dset["tau"])
print(dset["wp"])
print(dset["cldamt_ir"])
print(dset["cldamt_types"])
print(dset["snoice"])

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

#veo un primer grafico

cldamt_data=dset["cldamt"].mean("time", keep_attrs=True)

fig = plt.figure(figsize=[12,5])

ax = fig.add_subplot(111,projection=ccrs.PlateCarree(central_longitude=180))

cldamt_data.plot.contourf(ax=ax,
                   levels=np.arange(0, 100, 10),
                   extend='max',
                   transform=ccrs.PlateCarree(),
                   cbar_kwargs={'label': cldamt_data.units},
                   cmap='viridis_r')
ax.coastlines()

plt.show()


#%%
"""
Armo funcion que:
    Abre archivos nc 
    Selecciona region
    Selecciona variables
"""



