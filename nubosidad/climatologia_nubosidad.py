# climatologia_nubosidad.py

"""
Climatologia de nubosidad con datos de ISCCP-H HGM 

Los datos vienen dados como Promedios mensuales con una resolucion espacial de 1º, estan en formato netCDF CF complaint y se descargaron directamente o con THREDDS. Los detalles de los datos estan en este link: https://www.ncdc.noaa.gov/isccp/isccp-data-access/isccp-basic-data y la metadata aca: https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.ncdc:C00956 

Los promedios mensuales HGM surgen de hacer el promedio mensual de los datos ISCCP-H HGH que son los datos del promedio cada 3 horas para cada hora sinoptica. 

Los datos empiezan en el mes 07 del año 1983 y se extienden hasta el mes 06 del año 2017

Para descargar la informacion use la opcion de descarga directa con este comando en la terminal (fuente: https://stackoverflow.com/questions/6827459/can-i-use-wget-to-download-multiple-files-from-linux-terminal) :

    >wget -r -l1 -A.nc https://www.ncei.noaa.gov/data/international-satellite-cloud-climate-project-isccp-h-series-data/access/isccp-basic/hgm/

Uso este tutorial https://carpentrieslab.github.io/python-aos-lesson/02-visualisation/index.html para construir el codigo

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Pararse en la carpeta de Documentos/Doctorado/codigos/codigos2021/nubosidad 

Los datos estan en la carpeta Documentos/Doctorado/datos/nubosidad/ISCCP-H_HGM

"""

# %% Cargo librerias necesarias para el codigo

import geopandas
from cartopy.io import shapereader
from shapely.geometry.multipolygon import MultiPolygon
import cmocean
import pandas as pd
import xarray as xr
import numpy as np
from shapely.geometry import mapping, Polygon
import geopandas as gpd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import pymannkendall as mk

# %% Genero tres nuevas variables que van a representar el porcentaje de cldamt de nubes bajas medias y altas, siendo estos las sumas de cldamt_types : BAJAS:[0:6] MEDIAS:[6:12] ALTAS:[12:19]

for i in range(0, len(data_list)):

    cldamt_types = data_list[i]["cldamt_types"].values[0]

    bajas_array = xr.DataArray(data=np.reshape(sum(cldamt_types[:6]), (1, cldamt_types.shape[1], cldamt_types.shape[2])), dims=[
                               "time", "lat", "lon"]).assign_coords({"lon": data_list[i].lon.values, "lat": data_list[i].lat.values, "time": data_list[i].time.values})
    bajas_array.attrs["units"] = "percent"

    medias_array = xr.DataArray(data=np.reshape(sum(cldamt_types[6:12]), (1, cldamt_types.shape[1], cldamt_types.shape[2])), dims=[
                                "time", "lat", "lon"]).assign_coords({"lon": data_list[i].lon.values, "lat": data_list[i].lat.values, "time": data_list[i].time.values})
    medias_array.attrs["units"] = "percent"

    altas_array = xr.DataArray(data=np.reshape(sum(cldamt_types[12:18]), (1, cldamt_types.shape[1], cldamt_types.shape[2])), dims=[
                               "time", "lat", "lon"]).assign_coords({"lon": data_list[i].lon.values, "lat": data_list[i].lat.values, "time": data_list[i].time.values})
    altas_array.attrs["units"] = "percent"

    data_list[i] = data_list[i].assign(
        cldamt_bajas=bajas_array, cldamt_medias=medias_array, cldamt_altas=altas_array)

# elimino variables que ya no sirven
del(cldamt_types)
del(bajas_array)
del(medias_array)
del(altas_array)
del(i)

# %% Genero array3D donde cada capa es un tiempo: funcion. Esto se hace anual, por mes y por trimestre y va a ser la entrada a las funciones de calculo de estadisticos y tendencias


def array_3D(data_list, variable, mes_o_meses, clip):
    """
    Funcion que genera un array 3D con la informacion de una variable de los xarrays de la data_list

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, tomar lista en el rango de anios que se desean
    variable : str
        nombre variable
    mes_o_meses : None, str or tuple of str
        None: analisis anual
        str:numero de mes dos digitos
        tuple:tupla con meses como str
    clip: bool
        True: esta recortado por shape
        False: no esta recortado por shape
    Returns
    -------
    arr_3D : array
        array 3D con la informacion de cada elemento de data_list superpuesta.

    """
    if clip == False:
        if mes_o_meses == None:  # caso anual
            n = 1
            arr_3D = np.empty(data_list[0][variable].shape)
            for data in data_list:
                variable_data = data[variable].values
                arr_3D = np.concatenate([arr_3D, variable_data])
                n += 1
                arr_3D = np.reshape(
                    arr_3D, (n, data_list[0][variable].shape[1], data_list[0][variable].shape[2]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]

        elif type(mes_o_meses) == str:  # caso mensual
            n = 1
            arr_3D = np.empty(data_list[0][variable].shape)
            for data in data_list:
                if (str(data["time"].values[0])[5:7] == mes_o_meses):
                    variable_data = data[variable].values
                    arr_3D = np.concatenate([arr_3D, variable_data])
                    n = n+1
                    arr_3D = np.reshape(
                        arr_3D, (n, data_list[0][variable].shape[1], data_list[0][variable].shape[2]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]

        # caso trimestral. calculo la media por trimestre para cada anio VER
        elif type(mes_o_meses) == tuple:
            n = 1
            arr_3D = np.empty(data_list[0][variable].shape)
            for i in range(0, len(data_list)):
                if (str(data_list[i]["time"].values[0])[5:7] == mes_o_meses[0]):
                    variable_data1 = data_list[i][variable].values
                    variable_data2 = data_list[i+1][variable].values
                    variable_data3 = data_list[i+2][variable].values
                    variable_data_media = np.nanmean(
                        np.array([variable_data1, variable_data2, variable_data3]), axis=0)
                    arr_3D = np.concatenate([arr_3D, variable_data_media])
                    n += 1
                    arr_3D = np.reshape(
                        arr_3D, (n, data_list[0][variable].shape[1], data_list[0][variable].shape[2]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]
    elif clip == True:
        if mes_o_meses == None:  # caso anual
            n = 1
            arr_3D = np.empty((1,data_list[0].shape[0],data_list[0].shape[1]))
            for data in data_list:
                variable_data = data.values
                variable_data=variable_data.reshape(1,variable_data.shape[0],variable_data.shape[1])
                arr_3D = np.concatenate([arr_3D, variable_data])
                n += 1
                arr_3D = np.reshape(
                    arr_3D, (n, data_list[0].shape[0], data_list[0].shape[1]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]

        elif type(mes_o_meses) == str:  # caso mensual
            n = 1
            arr_3D = np.empty((1,data_list[0].shape[0],data_list[0].shape[1]))
            for data in data_list:
                if (data.name[-5:-3] == mes_o_meses):
                    variable_data = data.values
                    variable_data=variable_data.reshape(1,variable_data.shape[0],variable_data.shape[1])
                    arr_3D = np.concatenate([arr_3D, variable_data])
                    n = n+1
                    arr_3D = np.reshape(
                        arr_3D, (n, data_list[0].shape[0], data_list[0].shape[1]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]

        # caso trimestral. calculo la media por trimestre para cada anio VER
        elif type(mes_o_meses) == tuple:
            n = 1
            arr_3D = np.empty((1,data_list[0].shape[0],data_list[0].shape[1]))
            for i in range(0, len(data_list)):
                if (data_list[i].name[-5:-3] == mes_o_meses[0]):
                    variable_data1 = data_list[i].values
                    variable_data1=variable_data1.reshape(1,variable_data1.shape[0],variable_data1.shape[1])
                    variable_data2 = data_list[i+1].values
                    variable_data2=variable_data2.reshape(1,variable_data2.shape[0],variable_data2.shape[1])
                    variable_data3 = data_list[i+2].values
                    variable_data3=variable_data3.reshape(1,variable_data3.shape[0],variable_data3.shape[1])
                    variable_data_media = np.nanmean(
                        np.array([variable_data1, variable_data2, variable_data3]), axis=0)
                    arr_3D = np.concatenate([arr_3D, variable_data_media])
                    n += 1
                    arr_3D = np.reshape(
                        arr_3D, (n, data_list[0].shape[0], data_list[0].shape[1]))
            arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]

    return arr_3D
# %% Calculos estadisticos: funcion. Media y desvio climatologico anual, por mes y por estacion.
def media_desvio(data_list, variable, mes_o_meses, unidad, clip):
    """
    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, tomar lista en el rango de anios que se desean
    variable : str
        nombre variable
    mes_o_meses : None, str or tuple of str
        None: analisis anual
        str:numero de mes dos digitos
        tuple:tupla con meses como str
    clip: bool
        True: esta recortado por shape
        False: no esta recortado por shape
    Returns
    -------
    [0] xarray con media del mes seleccionado
    [1] xarray con desvio del mes seleccionado

    """
    arr_3D = array_3D(data_list, variable, mes_o_meses, clip)
    media = np.nanmean(arr_3D, axis=0)
    desvio = np.nanstd(arr_3D, axis=0)

    # cargo lats y lons para armar xarray
    lats = data_list[0].coords["lat"].values
    lons = data_list[0].coords["lon"].values
    coords = [("lat", lats), ("lon", lons)]

    # salida media
    if mes_o_meses == None:
        xarray_media = xr.DataArray(media, coords=coords)
        xarray_media.name = variable + " media anual"
        xarray_media.attrs['units'] = unidad
    elif type(mes_o_meses) == str:  # caso mensual
        xarray_media = xr.DataArray(media, coords=coords)
        xarray_media.name = variable + " media mensual " + mes_o_meses
        xarray_media.attrs['units'] = unidad
        xarray_media.attrs['mes'] = mes_o_meses
    elif type(mes_o_meses) == tuple:  # caso trimestral
        xarray_media = xr.DataArray(media, coords=coords)
        xarray_media.name = variable + " media trimestral " + \
            mes_o_meses[0]+"-"+mes_o_meses[1]+"-"+mes_o_meses[2]
        xarray_media.attrs['units'] = unidad
        xarray_media.attrs['meses'] = [
            mes_o_meses[0], mes_o_meses[1], mes_o_meses[2]]

    # salida desvio
    if mes_o_meses == None:
        xarray_desvio = xr.DataArray(desvio, coords=coords)
        xarray_desvio.name = variable + " desvio anual"
        xarray_desvio.attrs['units'] = unidad
    elif type(mes_o_meses) == str:  # caso mensual
        xarray_desvio = xr.DataArray(desvio, coords=coords)
        xarray_desvio.name = variable + " desvio mensual " + mes_o_meses
        xarray_desvio.attrs['units'] = unidad
        xarray_desvio.attrs['mes'] = mes_o_meses
    elif type(mes_o_meses) == tuple:  # caso trimestral
        xarray_desvio = xr.DataArray(desvio, coords=coords)
        xarray_desvio.name = variable + " desvio trimestral " + \
            mes_o_meses[0]+"-"+mes_o_meses[1]+"-"+mes_o_meses[2]
        xarray_desvio.attrs['units'] = unidad
        xarray_desvio.attrs['meses'] = [
            mes_o_meses[0], mes_o_meses[1], mes_o_meses[2]]

    return([xarray_media, xarray_desvio])
# %% Calculo tendencias con serie de datos: funcion. Calcula la tendencia y la testea con test de Mann Kendall


def tendencia(data):
    """

    Parameters
    ----------
    data : pd.DATAFRAME o np.array
        Data frame con serie de datos a calcular la tendencia. En la primer columna las fechas, en la segunda los dato

    Returns
    -------
    tendencia: NP.FLOAT
        valor de la tendencia dada entre un intervalo de tiempo. A esta salida si quiero tendencia decadal e ingrese datos de toda la serie mensual lo multiplico por 10*12 y si quiero por decada y vienen dados por un dato anual se multiplica por 10

    significativo: True False
        True si es significativo con un 95% de confianza con el test de Mann Kendall. False si no lo es. 
    """
    if type(data) == pd.core.frame.DataFrame:
        fechas = data.iloc[:, 0]
        fechas_list_1 = np.arange(0, len(fechas), 1)
        datos = data.iloc[:, 1]
        datos_array_1 = np.array(datos)
    elif type(data) == np.ndarray:
        fechas_list_1 = np.arange(0, len(data), 1)
        datos_array_1 = data

    datos_array = datos_array_1[np.isfinite(datos_array_1)]

    # si el vector es todo de nan entonces que no calcule la tendencia, y que devuelva nan
    if datos_array.size <= 1:
        tendencia = np.nan
        significativo = np.nan
    # si el vector tiene valores distintos de nan calcula tendencia y significancia
    elif datos_array.size > 1:
        fechas_list = fechas_list_1[np.isfinite(datos_array_1)]
        coef = np.polyfit(fechas_list, datos_array, 1)
        tendencia = coef[0]  # dada entre un intervalo de tiempo. A esta salida si quiero tendencia decadal e ingrese datos de toda la serie mensual lo multiplico por 10*12 y si quiero por decada y vienen dados por un dato anual se multiplica por 10
        z_test = mk.original_test(datos_array, alpha=0.05)[3]
        if abs(z_test) >= 1.96:
            significativo = True
        elif abs(z_test) < 1.96:
            significativo = False
    return(tendencia, significativo)
# %% Calculo tendencias para un campo: funcion. Calcula la tendencia decadal y la testea con test de Mann Kendall


def tendencia_campo(data_list, mes_o_meses, variable, unidad, pormesotrimestre, mes_trimestre, clip):
    """
    Calcula campo de tendencia

    Parameters
    ----------
    data_list : list
        lista en cada elemento un netcdf de un determinado mes y anio, tomar lista en el rango de anios que se desean
    mes_o_meses : None, str or tuple of str
        None: analisis anual
        str:numero de mes dos digitos
        tuple:tupla con meses como str
    variable : str
        nombre variable
    unidad : str
        unidad de la variable trabajada
    pormesotrimestre : bool
        False: cuando se quiere ver la tendencia anual
        True: cuando se quiere ver la tendencia por mes o por trimestre
    mes_trimestre : str
        Nombre de mes, de trimestre o "completo" para el analisis anual.
    clip : bool
        True: si la lista es de variable clipeada por el shape de corrientes
        False: si la lista es de variable no clipeada

    Returns
    -------
    None.

    """
    arr_3D = array_3D(data_list, variable, mes_o_meses,clip)

    if pormesotrimestre == False:
        n_lats = len(arr_3D[1, :, 1])
        n_lons = len(arr_3D[1, 1, :])
        tendencia_array = np.empty((2, n_lats, n_lons))
        for i in range(0, n_lats):
            for j in range(0, n_lons):
                data = arr_3D[:, i, j]
                tendencia_array[0, i, j] = tendencia(data)[0]*10*12
                if tendencia(data)[1] == True:
                    tendencia_array[1, i, j] = 1
                elif tendencia(data)[1] == False:
                    tendencia_array[1, i, j] = 0
    elif pormesotrimestre == True:
        n_lats = len(arr_3D[1, :, 1])
        n_lons = len(arr_3D[1, 1, :])
        tendencia_array = np.empty((2, n_lats, n_lons))
        for i in range(0, n_lats):
            for j in range(0, n_lons):
                data = arr_3D[:, i, j]
                tendencia_array[0, i, j] = tendencia(data)[0]*10
                if tendencia(data)[1] == True:
                    tendencia_array[1, i, j] = 1
                elif tendencia(data)[1] == False:
                    tendencia_array[1, i, j] = 0

    # cargo lats y lons para armar xarray
    if clip == False:
        lats = data_list[0].coords["lat"].values
        lons = data_list[0].coords["lon"].values
    if clip == True:
        lats = data_list[0].coords["y"].values
        lons = data_list[0].coords["x"].values
    coords = [("tendencia", np.array([1, 2])), ("lat", lats),
              ("lon", lons)]  # add first coord

    # salida tendencia
    xarray_media = xr.DataArray(tendencia_array, coords=coords)
    xarray_media.name = variable + " tendencia " + mes_trimestre
    xarray_media.attrs['units'] = unidad

    return(xarray_media)
    # return(tendencia_array)
#%% Recorto los xarray para el shape de Corrientes: funcion.
"""
Recorto los xarrays de la climatologia para el shape de Corrientes
Armo funcion que selecciona region de xarray con shape
https://gis.stackexchange.com/questions/382037/python-rioxarray-clip-masking-netcdf-data-with-a-polygon-returns-all-nan 
https://gis.stackexchange.com/questions/289775/masking-netcdf-data-using-shapefile-xarray-geopandas
"""


def clip_xarray_con_shape(coords, elemento_xarray_lista, climatologia_o_original, variable_original):
    """


    Parameters
    ----------
    coords : np.array
        Coordenadas del borde del shape que se quiere recortar. Por ejemplo: IGN=geopandas.read_file("/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp"); coords=np.array(IGN.geometry[19][0].exterior.coords)
    elemento_xarray_lista: list
        elemento de la lista con arrays que se quieren recortar.
    climatologia_o_original: str
        climatologia: para trabajar con xarrays de climatologia que tienen una unica variable
        original: para trabajar con xarrays originales
    variable_original: str
        variable que se esta trabajando, va a usarse solo en la funcion si se carga la lista con xarrays originales
    Returns
    -------
    List con primer elemento: xarray recortado, segundo elemento: shape que se recorto

    """
    coords_1=np.copy(coords)
    coords_1[:, 0] = coords_1[:, 0] +  360  # sumo 360 a las longitudes para que tenga la misma referencia que los xarrays
    newpoly = Polygon(coords_1)  # armo poligono
    newconus = gpd.GeoDataFrame(index=[0], crs="epsg:4326", geometry=[
                                newpoly])  # georeferencio el poligono

    # cargo xarray para recortar cambio el nombre de las dimensiones a x,y
    if (climatologia_o_original == "climatologia"):
        xds = elemento_xarray_lista
        xds_nombre = xds.name

    if (climatologia_o_original == "original"):
        xds = elemento_xarray_lista[variable_original]
        xds_nombre = xds.name+" "+str(xds["time"].values[0])[0:10]

    xds = xds.swap_dims({"lon": "x"})
    xds = xds.swap_dims({"lat": "y"})

    lats = xds["lat"][:].values
    lons = xds["lon"][:].values
    coords_2 = [("y", lats), ("x", lons)]

    if (climatologia_o_original == "climatologia"):
        xds = xr.DataArray(elemento_xarray_lista.values, coords=coords_2)
    if (climatologia_o_original == "original"):
        xds = xr.DataArray(
            elemento_xarray_lista[variable_original].values[0], coords=coords_2)

    xds = xds.rio.write_crs("epsg:4326", inplace=True)  # georeferencio

    # lo regrillo para que pase de ser 1ºx1º a 0.05ºx0.05º asi se recorta bien. Lo hago con el metodo de interpolacion lineal. el regrillado lo hago solo en un entorno a corrientes para ahorrar memoria
    # ynuevo=np.linspace(-89.5, 89.5, 3600) #si se cubriera todo el globo
    # xnuevo=np.linspace(0.5, 359.5, 7200) #si se cubriera todo el globo
    ynuevo = np.linspace(-31.5, -25.5, 60)
    xnuevo = np.linspace(280.5, 315.5, 350)
    # http://xarray.pydata.org/en/stable/generated/xarray.DataArray.interp.html
    xds = xds.interp(y=ynuevo, x=xnuevo, method="linear")

    clipped = xds.rio.clip(newconus.geometry.apply(mapping), newconus.crs, drop=False, invert=False)  # lo recorto
    clipped.name = xds_nombre
    # clipped.attrs['units']=xds.attrs['units']
    # clipped.attrs['mes']=xds.attrs['mes']
    # return(clipped,newconus)
    return clipped




# %% Calculo estadisticos: Corro funcion. Media y desvio anual, por mes y por estacion.


# defino lista con meses y con trimestres
meses = ["01", "02", "03", "04", "05", "06",
         "07", "08", "09", "10", "11", "12"]
trimestres = [("12", "01", "02"),("03", "04", "05"), ("06", "07", "08"),("09", "10", "11")]


# cldamt
# anuales
media_anual_cldamt = media_desvio(data_list[6:402], "cldamt", None, "%", False)[
    0]  # Enero 1984 - Diciembre 2016
desvio_anual_cldamt = media_desvio(data_list[6:402], "cldamt", None, "%", False)[
    1]  # Enero 1984 - Diciembre 2016

# medias y desvios mensuales
media_mensual_cldamt_list = []
desvio_mensual_cldamt_list = []
for mes in meses:
    media_mensual_cldamt_list.append(media_desvio(
        data_list[6:402], "cldamt", mes, "%", False)[0])  # Enero 1984 - Diciembre 2016
    desvio_mensual_cldamt_list.append(media_desvio(
        data_list[6:402], "cldamt", mes, "%", False)[1])  # Enero 1984 - Diciembre 2016

# medias y desvios trimestrales
media_trimestral_cldamt_list = []
desvio_trimestral_cldamt_list = []
for trimestre in trimestres:
    media_trimestral_cldamt_list.append(media_desvio(
        data_list[5:401], "cldamt", trimestre, "%", False)[0])  # Diciembre 1983 - Noviembre 2016
    desvio_trimestral_cldamt_list.append(media_desvio(
        data_list[5:401], "cldamt", trimestre, "%", False)[1])  # Diciembre 1983 - Noviembre 2016


# cldamt_bajas
# anuales
media_anual_cldamt_bajas = media_desvio(data_list[6:402], "cldamt_bajas", None, "%", False)[
    0]  # Enero 1984 - Diciembre 2016
desvio_anual_cldamt_bajas = media_desvio(
    data_list[6:402], "cldamt_bajas", None, "%", False)[1]  # Enero 1984 - Diciembre 2016

# medias y desvios mensuales
media_mensual_cldamt_bajas_list = []
desvio_mensual_cldamt_bajas_list = []
for mes in meses:
    media_mensual_cldamt_bajas_list.append(media_desvio(
        data_list[6:402], "cldamt_bajas", mes, "%",False)[0])  # Enero 1984 - Diciembre 2016
    desvio_mensual_cldamt_bajas_list.append(media_desvio(
        data_list[6:402], "cldamt_bajas", mes, "%",False)[1])  # Enero 1984 - Diciembre 2016

# medias y desvios trimestrales
media_trimestral_cldamt_bajas_list = []
desvio_trimestral_cldamt_bajas_list = []
for trimestre in trimestres:
    media_trimestral_cldamt_bajas_list.append(media_desvio(
        data_list[5:401], "cldamt_bajas", trimestre, "%",False)[0])  # Diciembre 1983 - Noviembre 2016
    desvio_trimestral_cldamt_bajas_list.append(media_desvio(
        data_list[5:401], "cldamt_bajas", trimestre, "%",False)[1])  # Diciembre 1983 - Noviembre 2016


# cldamt_medias
# anuales
media_anual_cldamt_medias = media_desvio(
    data_list[6:402], "cldamt_medias", None, "%",False)[0]  # Enero 1984 - Diciembre 2016
desvio_anual_cldamt_medias = media_desvio(
    data_list[6:402], "cldamt_medias", None, "%",False)[1]  # Enero 1984 - Diciembre 2016

# medias y desvios mensuales
media_mensual_cldamt_medias_list = []
desvio_mensual_cldamt_medias_list = []
for mes in meses:
    media_mensual_cldamt_medias_list.append(media_desvio(
        data_list[6:402], "cldamt_medias", mes, "%",False)[0])  # Enero 1984 - Diciembre 2016
    desvio_mensual_cldamt_medias_list.append(media_desvio(
        data_list[6:402], "cldamt_medias", mes, "%",False)[1])  # Enero 1984 - Diciembre 2016

# medias y desvios trimestrales
media_trimestral_cldamt_medias_list = []
desvio_trimestral_cldamt_medias_list = []
for trimestre in trimestres:
    media_trimestral_cldamt_medias_list.append(media_desvio(
        data_list[5:401], "cldamt_medias", trimestre, "%",False)[0])  # Diciembre 1983 - Noviembre 2016
    desvio_trimestral_cldamt_medias_list.append(media_desvio(
        data_list[5:401], "cldamt_medias", trimestre, "%",False)[1])  # Diciembre 1983 - Noviembre 2016


# cldamt_altas
# anuales
media_anual_cldamt_altas = media_desvio(data_list[6:402], "cldamt_altas", None, "%",False)[
    0]  # Enero 1984 - Diciembre 2016
desvio_anual_cldamt_altas = media_desvio(data_list[6:402], "cldamt_altas", None, "%",False)[
    1]  # Enero 1984 - Diciembre 2016

# medias y desvios mensuales
media_mensual_cldamt_altas_list = []
desvio_mensual_cldamt_altas_list = []
for mes in meses:
    media_mensual_cldamt_altas_list.append(media_desvio(
        data_list[6:402], "cldamt_altas", mes, "%",False)[0])  # Enero 1984 - Diciembre 2016
    desvio_mensual_cldamt_altas_list.append(media_desvio(
        data_list[6:402], "cldamt_altas", mes, "%",False)[1])  # Enero 1984 - Diciembre 2016

# medias y desvios trimestrales
media_trimestral_cldamt_altas_list = []
desvio_trimestral_cldamt_altas_list = []
for trimestre in trimestres:
    media_trimestral_cldamt_altas_list.append(media_desvio(
        data_list[5:401], "cldamt_altas", trimestre, "%",False)[0])  # Diciembre 1983 - Noviembre 2016
    desvio_trimestral_cldamt_altas_list.append(media_desvio(
        data_list[5:401], "cldamt_altas", trimestre, "%",False)[1])  # Diciembre 1983 - Noviembre 2016


# elimino variables que se generan que no van a servir luego
del(meses)
del(trimestres)
del(mes)
del(trimestre)

# %% Recorto los xarray para el shape de Corrientes: corro funcion. Para variable original y para los estadisticos media y desvio
"""
Corro esta funcion y armo lista con los estadisticos (media y desvio) recortados en la provincia de Corrientes mensual y trimestral y otra lista con todos los dias clippeados de cldamt
"""
# cargo shape con las provincia de argentina
# me muevo tres carpetas par atras y entro a datos
IGN = gpd.read_file("../../../datos/mapas/provincia/provincia.shp")

# cargo las coordenadas de corrientes
coords = np.array(IGN.geometry[19][0].exterior.coords)


# cldamt
# variable original
cldamt_list_clipped = []  # lista vacia
for data in data_list:
    cldamt_list_clipped.append(clip_xarray_con_shape(coords, data, "original", "cldamt"))
    
# media y desvio anual
media_anual_cldamt_clipped = clip_xarray_con_shape(
    coords, media_anual_cldamt, "climatologia", "cldamt")
desvio_anual_cldamt_clipped = clip_xarray_con_shape(
    coords, desvio_anual_cldamt, "climatologia", "cldamt")

# media y desvio mensual
media_mensual_cldamt_list_clipped = []
for data in media_mensual_cldamt_list:
    media_mensual_cldamt_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt"))
desvio_mensual_cldamt_list_clipped = []
for data in desvio_mensual_cldamt_list:
    desvio_mensual_cldamt_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt"))

# media y desvio trimestral
media_trimestral_cldamt_list_clipped = []
for data in media_trimestral_cldamt_list:
    media_trimestral_cldamt_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt"))
desvio_trimestral_cldamt_list_clipped = []
for data in desvio_trimestral_cldamt_list:
    desvio_trimestral_cldamt_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt"))


# cldamt_bajas
# variable original
cldamt_bajas_list_clipped = []  # lista vacia
for data in data_list:
    cldamt_bajas_list_clipped.append(clip_xarray_con_shape(
        coords, data, "original", "cldamt_bajas"))

# media y desvio anual
media_anual_cldamt_bajas_clipped = clip_xarray_con_shape(
    coords, media_anual_cldamt_bajas, "climatologia", "cldamt_bajas")
desvio_anual_cldamt_bajas_clipped = clip_xarray_con_shape(
    coords, desvio_anual_cldamt_bajas, "climatologia", "cldamt_bajas")

# media y desvio mensual
media_mensual_cldamt_bajas_list_clipped = []
for data in media_mensual_cldamt_bajas_list:
    media_mensual_cldamt_bajas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_bajas"))
desvio_mensual_cldamt_bajas_list_clipped = []
for data in desvio_mensual_cldamt_bajas_list:
    desvio_mensual_cldamt_bajas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_bajas"))

# media y desvio trimestral
media_trimestral_cldamt_bajas_list_clipped = []
for data in media_trimestral_cldamt_bajas_list:
    media_trimestral_cldamt_bajas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_bajas"))
desvio_trimestral_cldamt_bajas_list_clipped = []
for data in desvio_trimestral_cldamt_bajas_list:
    desvio_trimestral_cldamt_bajas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_bajas"))

# cldamt_medias
# variable original
cldamt_medias_list_clipped = []  # lista vacia
for data in data_list:
    cldamt_medias_list_clipped.append(clip_xarray_con_shape(
        coords, data, "original", "cldamt_medias"))

# media y desvio anual
media_anual_cldamt_medias_clipped = clip_xarray_con_shape(
    coords, media_anual_cldamt_medias, "climatologia", "cldamt_medias")
desvio_anual_cldamt_medias_clipped = clip_xarray_con_shape(
    coords, desvio_anual_cldamt_medias, "climatologia", "cldamt_medias")

# media y desvio mensual
media_mensual_cldamt_medias_list_clipped = []
for data in media_mensual_cldamt_medias_list:
    media_mensual_cldamt_medias_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_medias"))
desvio_mensual_cldamt_medias_list_clipped = []
for data in desvio_mensual_cldamt_medias_list:
    desvio_mensual_cldamt_medias_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_medias"))

# media y desvio trimestral
media_trimestral_cldamt_medias_list_clipped = []
for data in media_trimestral_cldamt_medias_list:
    media_trimestral_cldamt_medias_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_medias"))
desvio_trimestral_cldamt_medias_list_clipped = []
for data in desvio_trimestral_cldamt_medias_list:
    desvio_trimestral_cldamt_medias_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_medias"))

# cldamt_altas
# variable original
cldamt_altas_list_clipped = []  # lista vacia
for data in data_list:
    cldamt_altas_list_clipped.append(clip_xarray_con_shape(
        coords, data, "original", "cldamt_altas"))

# media y desvio anual
media_anual_cldamt_altas_clipped = clip_xarray_con_shape(
    coords, media_anual_cldamt_altas, "climatologia", "cldamt_altas")
desvio_anual_cldamt_altas_clipped = clip_xarray_con_shape(
    coords, desvio_anual_cldamt_altas, "climatologia", "cldamt_altas")

# media y desvio mensual
media_mensual_cldamt_altas_list_clipped = []
for data in media_mensual_cldamt_altas_list:
    media_mensual_cldamt_altas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_altas"))
desvio_mensual_cldamt_altas_list_clipped = []
for data in desvio_mensual_cldamt_altas_list:
    desvio_mensual_cldamt_altas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_altas"))

# media y desvio trimestral
media_trimestral_cldamt_altas_list_clipped = []
for data in media_trimestral_cldamt_altas_list:
    media_trimestral_cldamt_altas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_altas"))
desvio_trimestral_cldamt_altas_list_clipped = []
for data in desvio_trimestral_cldamt_altas_list:
    desvio_trimestral_cldamt_altas_list_clipped.append(
        clip_xarray_con_shape(coords, data, "climatologia", "cldamt_altas"))

# elimino variables que ya no me van a servir
del(IGN)
del(coords)
del(data)



#%% Calculo tendencias para un campo: corro funcion. Calcula la tendencia decadal y la testea con test de Mann Kendall

# defino lista con meses y con trimestres
meses = [("01","Enero"), ("02","Febrero"), ("03","Marzo"), ("04","Abril"), ("05","Mayo"), ("06","Junio"),("07","Julio"), ("08","Agosto"), ("09","Septiembre"), ("10","Octubre"), ("11","Noviembre"), ("12","Diciembre")]
trimestres = [(("12", "01", "02"),"DEF"),(("03", "04", "05"),"MAM"), (("06", "07", "08"),"JJA"),(("09", "10", "11"),"SON")]


# cldamt
# tendencia todo el anio, tendencia decadal y significancia
tendencia_anual_cldamt = tendencia_campo(data_list[6:402], None, "cldamt", "%/decada", False, "completo", False)  # Enero 1984 - Diciembre 2016
tendencia_anual_cldamt_clipped = tendencia_campo(cldamt_list_clipped[6:402], None, "cldamt", "%/decada", False, "completo", True)  # Enero 1984 - Diciembre 2016

# tendencia por mes, tendencia decadal y significancia
tendencia_mensual_cldamt_list = []
tendencia_mensual_cldamt_list_clipped = []
for mes in meses:
    tendencia_mensual_cldamt_list.append(tendencia_campo(data_list[6:402], mes[0], "cldamt", "%/decada", True ,mes[1], False))  # Enero 1984 - Diciembre 2016
    tendencia_mensual_cldamt_list_clipped.append(tendencia_campo(cldamt_list_clipped[6:402], mes[0], "cldamt", "%/decada", True ,mes[1], True))  # Enero 1984 - Diciembre 2016

# tendencia por trimestre, tendencia decadal y significancia
tendencia_trimestral_cldamt_list = []
tendencia_trimestral_cldamt_list_clipped = []
for trimestre in trimestres:
    tendencia_trimestral_cldamt_list.append(tendencia_campo(data_list[5:401], trimestre[0], "cldamt", "%/decada", True ,trimestre[1], False))  # Diciembre 1983 - Noviembre 2016
    tendencia_trimestral_cldamt_list_clipped.append(tendencia_campo(cldamt_list_clipped[5:401], trimestre[0], "cldamt", "%/decada", True ,trimestre[1], True))  # Diciembre 1983 - Noviembre 2016
   

# cldamt_bajas
# tendencia todo el anio, tendencia decadal y significancia
tendencia_anual_cldamt_bajas = tendencia_campo(data_list[6:402], None, "cldamt_bajas", "%/decada", False, "completo", False)  # Enero 1984 - Diciembre 2016
tendencia_anual_cldamt_bajas_clipped = tendencia_campo(cldamt_bajas_list_clipped[6:402], None, "cldamt_bajas", "%/decada", False, "completo", True)  # Enero 1984 - Diciembre 2016

# tendencia por mes, tendencia decadal y significancia
tendencia_mensual_cldamt_bajas_list = []
tendencia_mensual_cldamt_bajas_list_clipped = []
for mes in meses:
    tendencia_mensual_cldamt_bajas_list.append(tendencia_campo(data_list[6:402], mes[0], "cldamt_bajas", "%/decada", True ,mes[1], False))  # Enero 1984 - Diciembre 2016
    tendencia_mensual_cldamt_bajas_list_clipped.append(tendencia_campo(cldamt_bajas_list_clipped[6:402], mes[0], "cldamt_bajas", "%/decada", True ,mes[1], True))  # Enero 1984 - Diciembre 2016

# tendencia por trimestre, tendencia decadal y significancia
tendencia_trimestral_cldamt_bajas_list = []
tendencia_trimestral_cldamt_bajas_list_clipped = []
for trimestre in trimestres:
    tendencia_trimestral_cldamt_bajas_list.append(tendencia_campo(data_list[5:401], trimestre[0], "cldamt_bajas", "%/decada", True ,trimestre[1], False))  # Diciembre 1983 - Noviembre 2016
    tendencia_trimestral_cldamt_bajas_list_clipped.append(tendencia_campo(cldamt_bajas_list_clipped[5:401], trimestre[0], "cldamt_bajas", "%/decada", True ,trimestre[1], True))  # Diciembre 1983 - Noviembre 2016


# cldamt_medias
# tendencia todo el anio, tendencia decadal y significancia
tendencia_anual_cldamt_medias = tendencia_campo(data_list[6:402], None, "cldamt_medias", "%/decada", False, "completo", False)  # Enero 1984 - Diciembre 2016
tendencia_anual_cldamt_medias_clipped = tendencia_campo(cldamt_medias_list_clipped[6:402], None, "cldamt_medias", "%/decada", False, "completo", True)  # Enero 1984 - Diciembre 2016

# tendencia por mes, tendencia decadal y significancia
tendencia_mensual_cldamt_medias_list = []
tendencia_mensual_cldamt_medias_list_clipped = []
for mes in meses:
    tendencia_mensual_cldamt_medias_list.append(tendencia_campo(data_list[6:402], mes[0], "cldamt_medias", "%/decada", True ,mes[1], False))  # Enero 1984 - Diciembre 2016
    tendencia_mensual_cldamt_medias_list_clipped.append(tendencia_campo(cldamt_medias_list_clipped[6:402], mes[0], "cldamt_medias", "%/decada", True ,mes[1], True))  # Enero 1984 - Diciembre 2016

# tendencia por trimestre, tendencia decadal y significancia
tendencia_trimestral_cldamt_medias_list = []
tendencia_trimestral_cldamt_medias_list_clipped = []
for trimestre in trimestres:
    tendencia_trimestral_cldamt_medias_list.append(tendencia_campo(data_list[5:401], trimestre[0], "cldamt_medias", "%/decada", True ,trimestre[1], False))  # Diciembre 1983 - Noviembre 2016
    tendencia_trimestral_cldamt_medias_list_clipped.append(tendencia_campo(cldamt_medias_list_clipped[5:401], trimestre[0], "cldamt_medias", "%/decada", True ,trimestre[1], True))  # Diciembre 1983 - Noviembre 2016

# cldamt_altas
# tendencia todo el anio, tendencia decadal y significancia
tendencia_anual_cldamt_altas = tendencia_campo(data_list[6:402], None, "cldamt_altas", "%/decada", False, "completo", False)  # Enero 1984 - Diciembre 2016
tendencia_anual_cldamt_altas_clipped = tendencia_campo(cldamt_altas_list_clipped[6:402], None, "cldamt_altas", "%/decada", False, "completo", True)  # Enero 1984 - Diciembre 2016

# tendencia por mes, tendencia decadal y significancia
tendencia_mensual_cldamt_altas_list = []
tendencia_mensual_cldamt_altas_list_clipped = []
for mes in meses:
    tendencia_mensual_cldamt_altas_list.append(tendencia_campo(data_list[6:402], mes[0], "cldamt_altas", "%/decada", True ,mes[1], False))  # Enero 1984 - Diciembre 2016
    tendencia_mensual_cldamt_altas_list_clipped.append(tendencia_campo(cldamt_altas_list_clipped[6:402], mes[0], "cldamt_altas", "%/decada", True ,mes[1], True))  # Enero 1984 - Diciembre 2016

# tendencia por trimestre, tendencia decadal y significancia
tendencia_trimestral_cldamt_altas_list = []
tendencia_trimestral_cldamt_altas_list_clipped = []
for trimestre in trimestres:
    tendencia_trimestral_cldamt_altas_list.append(tendencia_campo(data_list[5:401], trimestre[0], "cldamt_altas", "%/decada", True ,trimestre[1], False))  # Diciembre 1983 - Noviembre 2016
    tendencia_trimestral_cldamt_altas_list_clipped.append(tendencia_campo(cldamt_altas_list_clipped[5:401], trimestre[0], "cldamt_altas", "%/decada", True ,trimestre[1], True))  # Diciembre 1983 - Noviembre 2016

#%%
#%% Grafico campos climatologia nubosidad: funcion. Falta corregir y unir con clip
"""
Defino funcion que grafica climatologia mensual de alguna variable definida previamente en una determinada region (region cuadrada)
"""


def grafico_campos_climatologia_nubosidad(paises, provincias, data_list_climatologia, indice_list, variable, climatologia_tipo, lat_min, lat_max, lon_min, lon_max, unidades_nombre, valor_minimo, valor_maximo, delta_valor, xticks_min, xticks_max, yticks_min, yticks_max, grid, region, ruta_salida, paleta_color, espacio_entre_lat_lon, orientacion):
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

    # limpio graficos
    plt.close()

    # selecciona variable
    variable_data = data_list_climatologia[indice_list]
    # variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    # cambio el nombre de la unidad
    variable_data.attrs["units"] = unidades_nombre

    # selecciono region
    lats = variable_data["lat"][:]
    lons = variable_data["lon"][:]
    lat_lims = [lat_min, lat_max]
    lon_lims = [360+lon_min, 360+lon_max]  # lean 360-64 (64 O) 360-31 (31 O)
    lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
    lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
    variable_data_subset = variable_data[lat_inds, lon_inds]

    # extraigo mes (climatologia mensual) o meses (climatologia trimestral)
    if (climatologia_tipo == "mensual"):
        meses = str(data_list_climatologia[indice_list].name)[-2] + \
            str(data_list_climatologia[indice_list].name)[-1]
        if (meses == "01"):
            mes = "Enero"
        if (meses == "02"):
            mes = "Febrero"
        if (meses == "03"):
            mes = "Marzo"
        if (meses == "04"):
            mes = "Abril"
        if (meses == "05"):
            mes = "Mayo"
        if (meses == "06"):
            mes = "Junio"
        if (meses == "07"):
            mes = "Julio"
        if (meses == "08"):
            mes = "Agosto"
        if (meses == "09"):
            mes = "Septiembre"
        if (meses == "10"):
            mes = "Octubre"
        if (meses == "11"):
            mes = "Noviembre"
        if (meses == "12"):
            mes = "Diciembre"

    if (climatologia_tipo == "trimestral"):
        meses = str(data_list_climatologia[indice_list].name)[-8:-1] + \
            str(data_list_climatologia[indice_list].name)[-1]
        if (meses == "03-04-05"):
            mes = "MAM"
        if (meses == "06-07-08"):
            mes = "JJA"
        if (meses == "09-10-11"):
            mes = "SON"
        if (meses == "12-01-02"):
            mes = "DEF"

    # ploteo
    if (orientacion == "H"):
        fig1 = plt.figure(figsize=[9, 5], dpi=200)  # horizontal region 1
    if (orientacion == "V"):
        fig1 = plt.figure(figsize=[7.5, 7.5], dpi=200)  # vertical sudamerica
    ax = fig1.add_subplot(
        111, projection=ccrs.PlateCarree(central_longitude=0))

    if (paleta_color == "rain"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.rain)

    if (paleta_color == "curl"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.curl_r)

    if (paleta_color == "matter"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.matter)
        #######ver estas lineas###############
    #facecolors = np.ma.array(variable_data_subset, mask=np.isnan(variable_data_subset))
    # ax.set_array(facecolors)
        #######ver estas lineas###############
    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.5', linewidth=0.7, alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.4', alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xticklabels(np.arange(xticks_min, xticks_max)
                       [::espacio_entre_lat_lon])
    plt.xticks(np.arange(xticks_min, xticks_max)[::espacio_entre_lat_lon])
    ax.set_xlabel("Longitud")

    ax.set_yticklabels(np.arange(yticks_min, yticks_max)
                       [::espacio_entre_lat_lon])
    plt.yticks(np.arange(yticks_min, yticks_max)[::espacio_entre_lat_lon])
    ax.set_ylabel("Latitud")

    if (grid == True):
        plt.grid(linestyle="--", alpha=0.3)

    plt.title(variable+" "+region+" "+mes)
    # plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes)
    plt.show()


#%% Grafico campos climatologia nubosidad clip: funcion. Falta corregir y unir con no clip
"""
Defino funcion que grafica climatologia mensual de alguna variable definida previamente en una determinada region (region con shape)
"""
# carga librerias necesarias


def grafico_campos_climatologia_nubosidad_clip(paises, provincias, data_list_climatologia, indice_list, variable, climatologia_tipo, lat_min, lat_max, lon_min, lon_max, unidades_nombre, valor_minimo, valor_maximo, delta_valor, xticks_min, xticks_max, yticks_min, yticks_max, grid, region, ruta_salida, paleta_color, espacio_entre_lat_lon, orientacion):
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
    # carga librerias necesarias
    #import matplotlib.pyplot as plt
    #import cartopy.crs as ccrs
    #import cmocean

    # limpio graficos
    plt.close()

    # selecciona variable
    variable_data = data_list_climatologia[indice_list]
    # variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    # variable_data.attrs["units"]=unidades_nombre #cambio el nombre de la unidad

    # selecciono region
    lats = variable_data["y"][:]
    lons = variable_data["x"][:]
    lat_lims = [lat_min, lat_max]
    lon_lims = [360+lon_min, 360+lon_max]  # lean 360-64 (64 O) 360-31 (31 O)
    lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
    lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
    variable_data_subset = variable_data[lat_inds, lon_inds]

    # extraigo mes (climatologia mensual) o meses (climatologia trimestral)
    if (climatologia_tipo == "mensual"):
        meses = str(data_list_climatologia[indice_list].name)[-2] + \
            str(data_list_climatologia[indice_list].name)[-1]
        if (meses == "01"):
            mes1 = "Enero"
        if (meses == "02"):
            mes1 = "Febrero"
        if (meses == "03"):
            mes1 = "Marzo"
        if (meses == "04"):
            mes1 = "Abril"
        if (meses == "05"):
            mes1 = "Mayo"
        if (meses == "06"):
            mes1 = "Junio"
        if (meses == "07"):
            mes1 = "Julio"
        if (meses == "08"):
            mes1 = "Agosto"
        if (meses == "09"):
            mes1 = "Septiembre"
        if (meses == "10"):
            mes1 = "Octubre"
        if (meses == "11"):
            mes1 = "Noviembre"
        if (meses == "12"):
            mes1 = "Diciembre"

    if (climatologia_tipo == "trimestral"):
        meses = str(data_list_climatologia[indice_list].name)[-8:-1] + \
            str(data_list_climatologia[indice_list].name)[-1]
        if (meses == "03-04-05"):
            mes1 = "MAM"
        if (meses == "06-07-08"):
            mes1 = "JJA"
        if (meses == "09-10-11"):
            mes1 = "SON"
        if (meses == "12-01-02"):
            mes1 = "DEF"

    # ploteo
    if (orientacion == "H"):
        fig1 = plt.figure(figsize=[9, 5], dpi=200)  # horizontal region 1
    if (orientacion == "V"):
        fig1 = plt.figure(figsize=[7.5, 7.5], dpi=200)  # vertical sudamerica
    ax = fig1.add_subplot(
        111, projection=ccrs.PlateCarree(central_longitude=0))

    if (paleta_color == "rain"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': unidades_nombre},
                                           cmap=cmocean.cm.rain)

    if (paleta_color == "curl"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': unidades_nombre},
                                           cmap=cmocean.cm.curl_r)

    if (paleta_color == "matter"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': unidades_nombre},
                                           cmap=cmocean.cm.matter)

    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.5', linewidth=0.7, alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.4', alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xticklabels(np.arange(xticks_min, xticks_max)
                       [::espacio_entre_lat_lon])
    plt.xticks(np.arange(xticks_min, xticks_max)[::espacio_entre_lat_lon])
    ax.set_xlabel("Longitud")

    ax.set_yticklabels(np.arange(yticks_min, yticks_max)
                       [::espacio_entre_lat_lon])
    plt.yticks(np.arange(yticks_min, yticks_max)[::espacio_entre_lat_lon])
    ax.set_ylabel("Latitud")

    if (grid == True):
        plt.grid(linestyle="--", alpha=0.3)

    #plt.title(variable+" "+region)
    plt.title(variable+" "+region+" "+mes1)
    # plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes1)
    plt.show()


#%% Grafico campos tendencias nubosidad: funcion. 
"""
Defino funcion que grafica climatologia de la tendencia definida previamente en una determinada region (region cuadrada)
# https://matplotlib.org/stable/gallery/images_contours_and_fields/irregulardatagrid.html
# https://stackoverflow.com/questions/36721977/overly-patches-which-represent-the-significants-points-over-contour-map
"""
# carga librerias necesarias
def grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_xarray, indice_list, variable, climatologia_tipo, lat_min, lat_max, lon_min, lon_max, unidades_nombre, valor_minimo, valor_maximo, delta_valor, xticks_min, xticks_max, yticks_min, yticks_max, grid, region, ruta_salida, paleta_color, espacio_entre_lat_lon, orientacion, nombre_periodo, color_fondo, periodo_titulo):
    """
    Parameters
    ----------
    paises : shapely.geometry.multipolygon.MultiPolygon
        shape con paises a graficar en mapa
    provincias : shapely.geometry.multipolygon.MultiPolygon
        shape con provincias a graficar en mapa
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
    # carga librerias necesarias
    #import matplotlib.pyplot as plt
    #import cartopy.crs as ccrs
    #import cmocean

    # limpio graficos
    plt.close()

    # selecciona variable
    variable_data = tendencia_xarray[0]
    variable_significancia = tendencia_xarray[1]
    # variable_data=data[variable].mean("time", keep_attrs=True) #selecciona variable y toma el unico valor para cada punto de grilla
    # variable_data.attrs["units"]=unidades_nombre #cambio el nombre de la unidad

    # selecciono region
    lats = variable_data["lat"][:]
    lons = variable_data["lon"][:]
    lat_lims = [lat_min, lat_max]
    lon_lims = [360+lon_min, 360+lon_max]  # lean 360-64 (64 O) 360-31 (31 O)
    lat_inds = np.where((lats.values > lat_lims[0]) & (
        lats.values < lat_lims[1]))[0]
    lon_inds = np.where((lons.values > lon_lims[0]) & (
        lons.values < lon_lims[1]))[0]
    variable_data_subset = variable_data[lat_inds, :][:, lon_inds]
    variable_significancia_subset = variable_significancia[lat_inds, :][:, lon_inds]

    # me quedo solo con los que No son significativos para marcar esas regiones con puntitos, si quiero marcar las significativas cambiar 0 por 1
    variable_significancia_subset = variable_significancia_subset.where(
        variable_significancia_subset.values == 1)

    # ploteo
    if (orientacion == "H"):
        fig1 = plt.figure(figsize=[9, 5], dpi=200)  # horizontal region 1
    if (orientacion == "V"):
        fig1 = plt.figure(figsize=[7.5, 7.5], dpi=200)  # vertical sudamerica
    ax = fig1.add_subplot(
        111, projection=ccrs.PlateCarree(central_longitude=0))

    if (paleta_color == "rain"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.rain)

    if (paleta_color == "curl"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.curl_r)

    if (paleta_color == "matter"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.matter)

    if (paleta_color == "balance"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='both',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.balance)

    # grafico puntos en los lugares donde no es significativo
    variable_significancia_subset.plot.contourf(ax=ax, levels=[0, 1], hatches=[
                                                '.'], alpha=0, add_colorbar=False)

    # seteo el fondo gris para que los nans aparezcan en este color
    if (color_fondo == "gris"):
        ax.set_facecolor('tab:gray')

    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.5', linewidth=0.7, alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.4', alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xlim(np.min(variable_data_subset["lon"][:]).values -
                360, np.max(variable_data_subset["lon"][:]).values-360)
    ax.set_xticklabels(np.arange(xticks_min, xticks_max)
                       [::espacio_entre_lat_lon])
    plt.xticks(np.arange(xticks_min, xticks_max)[::espacio_entre_lat_lon])
    ax.set_xlabel("Longitud")

    ax.set_ylim(np.min(variable_data_subset["lat"][:]).values, np.max(
        variable_data_subset["lat"][:]).values)
    ax.set_yticklabels(np.arange(yticks_min, yticks_max)
                       [::espacio_entre_lat_lon])
    plt.yticks(np.arange(yticks_min, yticks_max)[::espacio_entre_lat_lon])
    ax.set_ylabel("Latitud")

    if (grid == True):
        plt.grid(linestyle="--", alpha=0.3)

    plt.title(variable+" "+region+" "+periodo_titulo)
    # plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable +
                "_tendencia_"+region+"_"+nombre_periodo)
    plt.show()

#%% Grafico campos: cargo shapes paises y provincias. tendencias nubosidad: corro funcion. 

# cargo shape con paises

resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'
shpfilename = shapereader.natural_earth(
    resolution, category, name)  # cargo paises de natural_earth
# cargo el shapefile usando geopandas
df = geopandas.read_file(shpfilename)
# leo los paises que voy a necesitar
paises = MultiPolygon([df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == 'Brazil']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == 'Paraguay']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Uruguay']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Bolivia']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Chile']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Colombia"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Ecuador"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] ==
                              "Venezuela"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Guyana"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Suriname"]['geometry'].values[0],
                       df.loc[df['ADMIN'] == "Panama"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Costa Rica"]['geometry'].values[0][0]])  # los paso a multipolygon para poder graficarlos

# cargo shape con provincias de argentina con datos del IGN
# descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN = geopandas.read_file(
    "/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias = [None]*24
for i in range(0, 24):
    provincias[i] = IGN["geometry"][i]
# paso a multipolygon para poder ponerlo en mapa
provincias = MultiPolygon(provincias)

#elimino variables que no sirven
del(resolution)
del(category)
del(name)
del(shpfilename)
del(df)
del(IGN)

#%% Grafico campos tendencias nubosidad: corro funcion. 

#cargo lista con meses y estaciones
meses = ["Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio", "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"]
estaciones = ["DEF", "MAM", "JJA", "SON"]

# cldamt

# periodo completo
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt, 1, "cldamt tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 8, "V", "periodo_completo", "gris", "") # sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt, 1, "cldamt tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 4, "H", "periodo_completo", "gris", "") # region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt, 1, "cldamt tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 2, "H", "periodo_completo", "gris", "") # region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_clipped, 1, "cldamt tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "") # corrientes

plt.close(fig="all")

# por mes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 8, "V", meses[i], "gris", meses[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 4, "H", meses[i], "gris", meses[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 2, "H", meses[i], "gris", meses[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_list_clipped[i], 1, "cldamt tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i]) # corrientes

plt.close(fig="all")

# por estacion
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_list[i], 1, "cldamt tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8,9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_list_clipped[i], 1, "cldamt tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i]) # corrientes

plt.close(fig="all")

# cldamt_bajas

# periodo completo
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_bajas, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", "periodo_completo", "gris", "") # sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_bajas, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", "periodo_completo", "gris", "") # region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_bajas, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", "periodo_completo", "gris", "") # region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_bajas_clipped, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "") # corrientes

plt.close(fig="all")

# por mes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", meses[i], "gris", meses[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", meses[i], "gris", meses[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", meses[i], "gris", meses[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_bajas_list_clipped[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i]) # corrientes

plt.close(fig="all")

# por estacion
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_bajas_list[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8,9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_bajas_list_clipped[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i]) # corrientes

plt.close(fig="all")


# cldamt_medias

# periodo completo
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_medias, 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", "periodo_completo", "gris", "") # sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_medias, 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", "periodo_completo", "gris", "") # region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_medias, 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", "periodo_completo", "gris", "") # region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_medias_clipped, 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "") # corrientes

plt.close(fig="all")

# por mes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", meses[i], "gris", meses[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", meses[i], "gris", meses[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", meses[i], "gris", meses[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_medias_list_clipped[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i]) # corrientes

plt.close(fig="all")

# por estacion
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_medias_list[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8,9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_medias_list_clipped[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i]) # corrientes

plt.close(fig="all")


# cldamt_altas

# periodo completo
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_altas, 1, "cldamt altas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 8, "V", "periodo_completo", "gris", "") # sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_altas, 1, "cldamt altas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 4, "H", "periodo_completo", "gris", "") # region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_altas, 1, "cldamt altas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 2, "H", "periodo_completo", "gris", "") # region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_anual_cldamt_altas_clipped, 1, "cldamt altas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "") # corrientes

plt.close(fig="all")

# por mes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 8, "V", meses[i], "gris", meses[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 4, "H", meses[i], "gris", meses[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 2, "H", meses[i], "gris", meses[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_mensual_cldamt_altas_list_clipped[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i]) # corrientes

plt.close(fig="all")

# por estacion
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i]) # sudamerica
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -39, -16, -64, -31, "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i]) # region 1
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_altas_list[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -32, -22, -64, -53, "%", -8,9, 1, -63, -53, -31, -22, True, "Región 2", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i]) # region 2
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, tendencia_trimestral_cldamt_altas_list_clipped[i], 1, "cldamt altas tendencia decadal (1984-2016)", "", -31, -26, -60, -55, "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "../../../resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i]) # corrientes

plt.close(fig="all")
# %%

# cldamt bajas
# periodo completo
# ploteo para sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_COMPLETO_subset, TENDENCIA_BAJAS_COMPLETO, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8,
                                                 9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", "periodo_completo", "gris", "")
# region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_COMPLETO_subset, TENDENCIA_BAJAS_COMPLETO, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                 "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", "periodo_completo", "gris", "")
# region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_COMPLETO_subset, TENDENCIA_BAJAS_COMPLETO, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                 "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", "periodo_completo", "gris", "")
# corrientes
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_CLIPPED_COMPLETO, TENDENCIA_BAJAS_CLIP_COMPLETO, 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                 "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "")

plt.close(fig="all")

# por mes
meses = ["Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
         "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"]

# sudamerica
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_MESES_subset, TENDENCIA_BAJAS_MESES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8,
                                                     9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", meses[i], "gris", meses[i])

# region 1
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_MESES_subset, TENDENCIA_BAJAS_MESES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                     "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", meses[i], "gris", meses[i])

# region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_MESES_subset, TENDENCIA_BAJAS_MESES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                     "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", meses[i], "gris", meses[i])

# corrientes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_CLIPPED_MESES, TENDENCIA_BAJAS_CLIP_MESES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                     "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i])

# por estacion
estaciones = ["DEF", "MAM", "JJA", "SON"]

# sudamerica
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_ESTACIONES_subset, TENDENCIA_BAJAS_ESTACIONES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -60, 15, -90, -30, "%", -8,
                                                     9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i])

# region 1
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_ESTACIONES_subset, TENDENCIA_BAJAS_ESTACIONES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                     "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i])

# region 2
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_ESTACIONES_subset, TENDENCIA_BAJAS_ESTACIONES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                     "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i])

# corrientes
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_BAJAS_3D_CLIPPED_ESTACIONES, TENDENCIA_BAJAS_CLIP_ESTACIONES[i], 1, "cldamt bajas tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                     "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i])

# %%

# cldamt medias
# periodo completo
# ploteo para sudamerica
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_COMPLETO_subset, TENDENCIA_MEDIAS_COMPLETO, 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30,
                                                 "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", "periodo_completo", "gris", "")
# region 1
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_COMPLETO_subset, TENDENCIA_MEDIAS_COMPLETO, 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                 "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", "periodo_completo", "gris", "")
# region 2
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_COMPLETO_subset, TENDENCIA_MEDIAS_COMPLETO, 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                 "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", "periodo_completo", "gris", "")
# corrientes
grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_CLIPPED_COMPLETO, TENDENCIA_MEDIAS_CLIP_COMPLETO, 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                 "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", "periodo_completo", "blanco", "")

plt.close(fig="all")

# por mes
meses = ["Enero", "Febrero", "Marzo", "Abril", "Mayo", "Junio",
         "Julio", "Agosto", "Septiembre", "Octubre", "Noviembre", "Diciembre"]

# sudamerica
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_MESES_subset, TENDENCIA_MEDIAS_MESES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30,
                                                     "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", meses[i], "gris", meses[i])

# region 1 HACER EN 5,12
for i in range(5, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_MESES_subset, TENDENCIA_MEDIAS_MESES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                     "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", meses[i], "gris", meses[i])

# region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_MESES_subset, TENDENCIA_MEDIAS_MESES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                     "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", meses[i], "gris", meses[i])

# corrientes
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_CLIPPED_MESES, TENDENCIA_MEDIAS_CLIP_MESES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                     "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", meses[i], "blanco", meses[i])

# por estacion
estaciones = ["DEF", "MAM", "JJA", "SON"]

# sudamerica
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_ESTACIONES_subset, TENDENCIA_MEDIAS_ESTACIONES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -60, 15, -90, -30,
                                                     "%", -8, 9, 1, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 8, "V", estaciones[i], "gris", estaciones[i])

# region 1
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_ESTACIONES_subset, TENDENCIA_MEDIAS_ESTACIONES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -39, -16, -64, -31,
                                                     "%", -8, 9, 1, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 4, "H", estaciones[i], "gris", estaciones[i])

# region 2
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_ESTACIONES_subset, TENDENCIA_MEDIAS_ESTACIONES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -32, -22, -64, -53,
                                                     "%", -8, 9, 1, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 2, "H", estaciones[i], "gris", estaciones[i])

# corrientes
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_tendencias(paises, provincias, ARR_MEDIAS_3D_CLIPPED_ESTACIONES, TENDENCIA_MEDIAS_CLIP_ESTACIONES[i], 1, "cldamt medias tendencia decadal (1984-2016)", "", -31, -26, -60, -55,
                                                     "%", -8, 9, 1, -59, -55, -30, -26, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "balance", 1, "H", estaciones[i], "blanco", estaciones[i])

"""
Extender todo a altas
"""
# %%
"""
Grafico climatologias
"""


# cargo shape con paises

resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'
shpfilename = shapereader.natural_earth(
    resolution, category, name)  # cargo paises de natural_earth
# cargo el shapefile usando geopandas
df = geopandas.read_file(shpfilename)
# leo los paises que voy a necesitar
paises = MultiPolygon([df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == 'Brazil']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == 'Paraguay']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Uruguay']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Bolivia']['geometry'].values[0],
                       df.loc[df['ADMIN'] == 'Chile']['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Colombia"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Ecuador"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] ==
                              "Venezuela"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Guyana"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Suriname"]['geometry'].values[0],
                       df.loc[df['ADMIN'] == "Panama"]['geometry'].values[0][0],
                       df.loc[df['ADMIN'] == "Costa Rica"]['geometry'].values[0][0]])  # los paso a multipolygon para poder graficarlos

# cargo shape con provincias de argentina con datos del IGN
# descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN = geopandas.read_file(
    "/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias = [None]*24
for i in range(0, 24):
    provincias[i] = IGN["geometry"][i]
# paso a multipolygon para poder ponerlo en mapa
provincias = MultiPolygon(provincias)

# %%
# ploteo para sudamerica
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_list, i, "cldamt media mensual (1984-2016)", "mensual", -60, 15, -90, -30, "%",
                                          0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_list, i, "cldamt desvío estándar mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 26, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 8, "V")

# otra escala la misma que para cldamt medias y altas
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_list, i, "cldamt desvío estándar mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica ", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 8, "V")


plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_list, i, "cldamt media trimestral (1984-2016)", "trimestral", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_list, i, "cldamt desvío estándar trimestral (1984-2016)", "trimestral", -60, 15, -90, -30,
                                          "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica ", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 8, "V")

plt.close(fig="all")

# ploteo para region 1
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_list, i, "cldamt media mensual (1984-2016)", "mensual", -39, -16, -64, -31, "%",
                                          0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_list, i, "cldamt desvío estándar mensual (1984-2016)", "mensual", -39, -16, -64, -31,
                                          "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_list, i, "cldamt media trimestral (1984-2016)", "trimestral", -39, -16, -64, -31,
                                          "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_list, i, "cldamt desvío estándar trimestral (1984-2016)", "trimestral", -39, -16, -64, -31,
                                          "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 4, "H")

# ploteo para region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_list, i, "cldamt media mensual (1984-2016)", "mensual", -32, -22, -64, -53, "%",
                                          0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_list, i, "cldamt desvío estándar mensual (1984-2016)", "mensual", -32, -22, -64, -53,
                                          "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_list, i, "cldamt media trimestral (1984-2016)", "trimestral", -32, -22, -64, -53,
                                          "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_list, i, "cldamt desvío estándar trimestral (1984-2016)", "trimestral", -32, -22, -64, -53,
                                          "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 2, "H")

# PLOTEO PARA CORRIENTES
# plt.close(fig="all")

# paleta1
# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_list_clipped,i,"cldamt media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,101,5,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",1,"H")

plt.close(fig="all")

# paleta2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_mensual_cldamt_list_clipped, i, "cldamt media mensual (1984-2016)", "mensual", -31, -26, -60, -55,
                                               "%", 40, 70, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_list_clipped,i,"cldamt desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,26,2,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",1,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_mensual_cldamt_list_clipped, i, "cldamt desvío estándar mensual (1984-2016)", "mensual", -31, -26, -
                                               60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 1, "H")

# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_list,i,"cldamt media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","rain",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_trimestral_cldamt_list_clipped, i, "cldamt media trimestral (1984-2016)", "trimestral", -31, -26, -
                                               60, -55, "%", 40, 70, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_list,i,"cldamt desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia","matter",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_trimestral_cldamt_list_clipped, i, "cldamt desvío estándar trimestral (1984-2016)", "trimestral", -31, -
                                               26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_climatologia", "matter", 1, "H")


# %%
# extiendo a cldamt bajas medias y altas
# %%
# bajas
# ploteo para sudamerica
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_bajas_list, i, "cldamt nubes bajas media mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar mensual (1984-2016)", "mensual", -60, 15, -
                                          90, -30, "%", 0, 26, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 8, "V")

# cambio escala
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica ", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 8, "V")


plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas media trimestral (1984-2016)", "trimestral", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar trimestral (1984-2016)", "trimestral", -60, 15, -
                                          90, -30, "%", 0, 26, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 8, "V")

# cambio escala
for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar trimestral (1984-2016)", "trimestral", -60, 15, -
                                          90, -30, "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica ", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 8, "V")

plt.close(fig="all")

# region 1
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_bajas_list, i, "cldamt nubes bajas media mensual (1984-2016)", "mensual", -39, -16, -64, -31,
                                          "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar mensual (1984-2016)", "mensual", -39, -16, -
                                          64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas media trimestral (1984-2016)", "trimestral", -39, -16, -
                                          64, -31, "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar trimestral (1984-2016)", "trimestral", -39, -
                                          16, -64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 4, "H")

# ploteo para region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_bajas_list, i, "cldamt nubes bajas media mensual (1984-2016)", "mensual", -32, -22, -64, -53,
                                          "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar mensual (1984-2016)", "mensual", -32, -22, -
                                          64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas media trimestral (1984-2016)", "trimestral", -32, -22, -
                                          64, -53, "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_bajas_list, i, "cldamt nubes bajas desvío estándar trimestral (1984-2016)", "trimestral", -32, -
                                          22, -64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 2, "H")

# PLOTEO PARA CORRIENTES
# plt.close(fig="all")

# paleta1
# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_bajas_list_clipped,i,"cldamt nubes bajas media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,101,5,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","rain",1,"H")

plt.close(fig="all")

# paleta2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_mensual_cldamt_bajas_list_clipped, i, "cldamt nubes bajas media mensual (1984-2016)", "mensual", -31, -26, -
                                               60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_bajas_list_clipped,i,"cldamt nubes bajas desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,26,2,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","matter",1,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_mensual_cldamt_bajas_list_clipped, i, "cldamt nubes bajas desvío estándar mensual (1984-2016)", "mensual", -31, -
                                               26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 1, "H")

# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_bajas_list,i,"cldamt nubes bajas media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","rain",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_trimestral_cldamt_bajas_list_clipped, i, "cldamt nubes bajas media trimestral (1984-2016)", "trimestral", -31, -
                                               26, -60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_bajas_list,i,"cldamt nubes bajas desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia","matter",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_trimestral_cldamt_bajas_list_clipped, i, "cldamt nubes bajas desvío estándar trimestral (1984-2016)", "trimestral", -
                                               31, -26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_climatologia", "matter", 1, "H")

# %%
###############################################################################################################################
# medias
# SUDAMERICA
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_medias_list, i, "cldamt nubes medias media mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_medias_list, i, "cldamt nubes medias desvío estándar mensual (1984-2016)", "mensual", -60, 15, -
                                          90, -30, "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_medias_list, i, "cldamt nubes medias media trimestral (1984-2016)", "trimestral", -60, 15, -
                                          90, -30, "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_medias_list, i, "cldamt nubes medias desvío estándar trimestral (1984-2016)", "trimestral", -60,
                                          15, -90, -30, "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 8, "V")

plt.close(fig="all")


# region 1
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_medias_list, i, "cldamt nubes medias media mensual (1984-2016)", "mensual", -39, -16, -64, -31,
                                          "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_medias_list, i, "cldamt nubes medias desvío estándar mensual (1984-2016)", "mensual", -39, -16, -
                                          64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_medias_list, i, "cldamt nubes medias media trimestral (1984-2016)", "trimestral", -39, -16, -
                                          64, -31, "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_medias_list, i, "cldamt nubes medias desvío estándar trimestral (1984-2016)", "trimestral", -39, -
                                          16, -64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 4, "H")

# ploteo para region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_medias_list, i, "cldamt nubes medias media mensual (1984-2016)", "mensual", -32, -22, -64, -53,
                                          "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_medias_list, i, "cldamt nubes medias desvío estándar mensual (1984-2016)", "mensual", -32, -22, -
                                          64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_medias_list, i, "cldamt nubes medias media trimestral (1984-2016)", "trimestral", -32, -22, -
                                          64, -53, "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_medias_list, i, "cldamt nubes medias desvío estándar trimestral (1984-2016)", "trimestral", -32, -
                                          22, -64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 2, "H")

# PLOTEO PARA CORRIENTES
# plt.close(fig="all")

# paleta1
# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_medias_list_clipped,i,"cldamt nubes medias media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,101,5,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","rain",1,"H")

plt.close(fig="all")

# paleta2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_mensual_cldamt_medias_list_clipped, i, "cldamt nubes medias media mensual (1984-2016)", "mensual", -31, -26, -
                                               60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_medias_list_clipped,i,"cldamt nubes medias desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,26,2,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","matter",1,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_mensual_cldamt_medias_list_clipped, i, "cldamt nubes medias desvío estándar mensual (1984-2016)", "mensual", -31, -
                                               26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 1, "H")

# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_medias_list,i,"cldamt nubes medias media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","rain",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_trimestral_cldamt_medias_list_clipped, i, "cldamt nubes medias media trimestral (1984-2016)", "trimestral", -31, -
                                               26, -60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_medias_list,i,"cldamt nubes medias desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia","matter",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_trimestral_cldamt_medias_list_clipped, i, "cldamt nubes medias desvío estándar trimestral (1984-2016)", "trimestral", -
                                               31, -26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_climatologia", "matter", 1, "H")

# %%
###############################################################################################################################
# altas
# SUDAMERICA

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_altas_list, i, "cldamt nubes altas media mensual (1984-2016)", "mensual", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 8, "V")
# ver que quedan espacios en blanco.

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_altas_list, i, "cldamt nubes altas desvío estándar mensual (1984-2016)", "mensual", -60, 15, -
                                          90, -30, "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_altas_list, i, "cldamt nubes altas media trimestral (1984-2016)", "trimestral", -60, 15, -90, -30,
                                          "%", 0, 101, 5, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 8, "V")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_altas_list, i, "cldamt nubes altas desvío estándar trimestral (1984-2016)", "trimestral", -60, 15, -
                                          90, -30, "%", 0, 40, 2, -85, -30, -55, 15, True, "Sudamérica", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 8, "V")

plt.close(fig="all")


# region 1
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_altas_list, i, "cldamt nubes altas media mensual (1984-2016)", "mensual", -39, -16, -64, -31,
                                          "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_altas_list, i, "cldamt nubes altas desvío estándar mensual (1984-2016)", "mensual", -39, -16, -
                                          64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_altas_list, i, "cldamt nubes altas media trimestral (1984-2016)", "trimestral", -39, -16, -
                                          64, -31, "%", 0, 101, 5, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 4, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_altas_list, i, "cldamt nubes altas desvío estándar trimestral (1984-2016)", "trimestral", -39, -
                                          16, -64, -31, "%", 0, 26, 2, -60, -31, -35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 4, "H")

# ploteo para region 2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_mensual_cldamt_altas_list, i, "cldamt nubes altas media mensual (1984-2016)", "mensual", -32, -22, -64, -53,
                                          "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_mensual_cldamt_altas_list, i, "cldamt nubes altas desvío estándar mensual (1984-2016)", "mensual", -32, -22, -
                                          64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, media_trimestral_cldamt_altas_list, i, "cldamt nubes altas media trimestral (1984-2016)", "trimestral", -32, -22, -
                                          64, -53, "%", 0, 101, 5, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 2, "H")

plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad(paises, provincias, desvio_trimestral_cldamt_altas_list, i, "cldamt nubes altas desvío estándar trimestral (1984-2016)", "trimestral", -32, -
                                          22, -64, -53, "%", 0, 26, 2, -63, -53, -31, -22, True, "Región 2", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 2, "H")

# PLOTEO PARA CORRIENTES
# plt.close(fig="all")

# paleta1
# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,media_mensual_cldamt_altas_list_clipped,i,"cldamt nubes altas media mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,101,5,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","rain",1,"H")

plt.close(fig="all")

# paleta2
for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_mensual_cldamt_altas_list_clipped, i, "cldamt nubes altas media mensual (1984-2016)", "mensual", -31, -26, -
                                               60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,12):
#    grafico_campos_climatologia_nubosidad_clip(paises,provincias,desvio_mensual_cldamt_altas_list_clipped,i,"cldamt nubes altas desvío estándar mensual (1984-2016)","mensual",-31,-26,-60,-55,"%",0,26,2,-59,-55,-30,-27,True,"Corrientes","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","matter",1,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 12):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_mensual_cldamt_altas_list_clipped, i, "cldamt nubes altas desvío estándar mensual (1984-2016)", "mensual", -31, -
                                               26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 1, "H")

# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,media_trimestral_cldamt_altas_list,i,"cldamt nubes altas media trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,101,5,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","rain",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, media_trimestral_cldamt_altas_list_clipped, i, "cldamt nubes altas media trimestral (1984-2016)", "trimestral", -31, -
                                               26, -60, -55, "%", 0, 40, 2, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "rain", 1, "H")


# plt.close(fig="all")

# for i in range(0,4):
#    grafico_campos_climatologia_nubosidad(paises,provincias,desvio_trimestral_cldamt_altas_list,i,"cldamt nubes altas desvío estándar trimestral (1984-2016)","trimestral",-39,-16,-64,-31,"%",0,26,2,-60,-31,-35,-18,True,"Región 1","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia","matter",4,"H")

# paleta2
plt.close(fig="all")

for i in range(0, 4):
    grafico_campos_climatologia_nubosidad_clip(paises, provincias, desvio_trimestral_cldamt_altas_list_clipped, i, "cldamt nubes altas desvío estándar trimestral (1984-2016)", "trimestral", -
                                               31, -26, -60, -55, "%", 0, 14, 1, -59, -55, -30, -27, True, "Corrientes", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_climatologia", "matter", 1, "H")

# %%

"""
calculo la media mensual de una determinada variable en determinada region 
grafico una serie temporal de todo el periodo, y por mes
"""
# defino funcion que calcula media mensual de una determinada variable en determinada region


def media_espacial(data_list, indice_list, variable, lat_min, lat_max, lon_min, lon_max, clipped):
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

    data = data_list[indice_list]
    if (clipped == False):
        # selecciona variable y toma el unico valor para cada punto de grilla
        variable_data = data[variable].mean("time", keep_attrs=True)
        # selecciono region
        lats = variable_data["lat"][:]
        lons = variable_data["lon"][:]
        lat_lims = [lat_min, lat_max]
        # lean 360-64 (64 O) 360-31 (31 O)
        lon_lims = [360+lon_min, 360+lon_max]
        lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
        lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
        variable_data_subset = variable_data[lat_inds, lon_inds]

    if (clipped == True):
        variable_data = data
        # selecciono region
        lats = variable_data["y"][:]
        lons = variable_data["x"][:]
        lat_lims = [lat_min, lat_max]
        # lean 360-64 (64 O) 360-31 (31 O)
        lon_lims = [360+lon_min, 360+lon_max]
        lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
        lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
        variable_data_subset = variable_data[lat_inds, lon_inds]

    # calculo media espacial del subset de la variable data
    media_espacial = np.nanmean(variable_data_subset.values, axis=(0, 1))
    return(media_espacial)

# armo funcion que devuelve data frame donde primera columna sea la fecha y la segunda columna sea la media espacial


def media_espacial_df(data_list, variable, lat_min, lat_max, lon_min, lon_max, clipped):
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
    media_espacial_df = pd.DataFrame(
        columns=["fecha", "Media_espacial_"+variable])
    if (clipped == False):
        for i in range(0, len(data_list)):
            # extraigo mes
            mes = str(data_list[i]["time"].values[0])[5:7]
            # extraigo anio
            anio = str(data_list[i]["time"].values[0])[0:4]
            # ¢alculo media
            media_espacial_i = media_espacial(
                data_list, i, variable, lat_min, lat_max, lon_min, lon_max, False)
            media_espacial_df = media_espacial_df.append(
                {"fecha": mes+"-"+anio, "Media_espacial_"+variable: media_espacial_i}, ignore_index=True)

    if (clipped == True):
        for i in range(0, len(data_list)):
            # extraigo mes
            mes = str(data_list[i].name)[-5:-3]
            # extraigo anio
            anio = str(data_list[i].name)[-10:-6]
            # ¢alculo media
            media_espacial_i = media_espacial(
                data_list, i, variable, lat_min, lat_max, lon_min, lon_max, True)
            media_espacial_df = media_espacial_df.append(
                {"fecha": mes+"-"+anio, "Media_espacial_"+variable: media_espacial_i}, ignore_index=True)

    media_espacial_df["fecha"] = pd.to_datetime(media_espacial_df["fecha"])
    return(media_espacial_df)


# %% corro para cldamt

cldamt_media_espacial_df_sudamerica = media_espacial_df(
    data_list, "cldamt", -60, 15, -90, -30, False)
cldamt_media_espacial_df_region1 = media_espacial_df(
    data_list, "cldamt", -39, -16, -64, -31, False)
cldamt_media_espacial_df_region2 = media_espacial_df(
    data_list, "cldamt", -32, -22, -64, -53, False)
cldamt_media_espacial_df_corrientes = media_espacial_df(
    cldamt_list_clipped, "cldamt", -32, -22, -64, -53, True)

cldamt_media_espacial_df_sudamerica.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_sudamerica.csv")
cldamt_media_espacial_df_region1.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_region1.csv")
cldamt_media_espacial_df_region2.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_region2.csv")
cldamt_media_espacial_df_corrientes.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_corrientes.csv")

# %% corro para cldamt bajas
cldamt_bajas_media_espacial_df_sudamerica = media_espacial_df(
    data_list, "cldamt_bajas", -60, 15, -90, -30, False)
cldamt_bajas_media_espacial_df_region1 = media_espacial_df(
    data_list, "cldamt_bajas", -39, -16, -64, -31, False)
cldamt_bajas_media_espacial_df_region2 = media_espacial_df(
    data_list, "cldamt_bajas", -32, -22, -64, -53, False)
cldamt_bajas_media_espacial_df_corrientes = media_espacial_df(
    cldamt_bajas_list_clipped, "cldamt_bajas", -32, -22, -64, -53, True)

cldamt_bajas_media_espacial_df_sudamerica.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_sudamerica.csv")
cldamt_bajas_media_espacial_df_region1.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_region1.csv")
cldamt_bajas_media_espacial_df_region2.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_region2.csv")
cldamt_bajas_media_espacial_df_corrientes.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_corrientes.csv")

# %% corro para cldamt medias
cldamt_medias_media_espacial_df_sudamerica = media_espacial_df(
    data_list, "cldamt_medias", -60, 15, -90, -30, False)
cldamt_medias_media_espacial_df_region1 = media_espacial_df(
    data_list, "cldamt_medias", -39, -16, -64, -31, False)
cldamt_medias_media_espacial_df_region2 = media_espacial_df(
    data_list, "cldamt_medias", -32, -22, -64, -53, False)
cldamt_medias_media_espacial_df_corrientes = media_espacial_df(
    cldamt_medias_list_clipped, "cldamt_medias", -32, -22, -64, -53, True)

cldamt_medias_media_espacial_df_sudamerica.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_sudamerica.csv")
cldamt_medias_media_espacial_df_region1.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_region1.csv")
cldamt_medias_media_espacial_df_region2.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_region2.csv")
cldamt_medias_media_espacial_df_corrientes.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_corrientes.csv")

# %% corro para cldamt altas
cldamt_altas_media_espacial_df_sudamerica = media_espacial_df(
    data_list, "cldamt_altas", -60, 15, -90, -30, False)
cldamt_altas_media_espacial_df_region1 = media_espacial_df(
    data_list, "cldamt_altas", -39, -16, -64, -31, False)
cldamt_altas_media_espacial_df_region2 = media_espacial_df(
    data_list, "cldamt_altas", -32, -22, -64, -53, False)
cldamt_altas_media_espacial_df_corrientes = media_espacial_df(
    cldamt_altas_list_clipped, "cldamt_altas", -32, -22, -64, -53, True)

cldamt_altas_media_espacial_df_sudamerica.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_sudamerica.csv")
cldamt_altas_media_espacial_df_region1.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_region1.csv")
cldamt_altas_media_espacial_df_region2.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_region2.csv")
cldamt_altas_media_espacial_df_corrientes.to_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_corrientes.csv")

# %% ahora que estan guardados los df, los cargo directamente.

# cldamt
cldamt_media_espacial_df_sudamerica = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_sudamerica.csv", index_col=0)
cldamt_media_espacial_df_sudamerica["fecha"] = pd.to_datetime(
    cldamt_media_espacial_df_sudamerica["fecha"])
cldamt_media_espacial_df_region1 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_region1.csv", index_col=0)
cldamt_media_espacial_df_region1["fecha"] = pd.to_datetime(
    cldamt_media_espacial_df_region1["fecha"])
cldamt_media_espacial_df_region2 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_region2.csv", index_col=0)
cldamt_media_espacial_df_region2["fecha"] = pd.to_datetime(
    cldamt_media_espacial_df_region2["fecha"])
cldamt_media_espacial_df_corrientes = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_corrientes.csv", index_col=0)
cldamt_media_espacial_df_corrientes["fecha"] = pd.to_datetime(
    cldamt_media_espacial_df_corrientes["fecha"])

# cldamt_bajas
cldamt_bajas_media_espacial_df_sudamerica = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_sudamerica.csv", index_col=0)
cldamt_bajas_media_espacial_df_sudamerica["fecha"] = pd.to_datetime(
    cldamt_bajas_media_espacial_df_sudamerica["fecha"])
cldamt_bajas_media_espacial_df_region1 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_region1.csv", index_col=0)
cldamt_bajas_media_espacial_df_region1["fecha"] = pd.to_datetime(
    cldamt_bajas_media_espacial_df_region1["fecha"])
cldamt_bajas_media_espacial_df_region2 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_region2.csv", index_col=0)
cldamt_bajas_media_espacial_df_region2["fecha"] = pd.to_datetime(
    cldamt_bajas_media_espacial_df_region2["fecha"])
cldamt_bajas_media_espacial_df_corrientes = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_bajas_media_espacial_df_corrientes.csv", index_col=0)
cldamt_bajas_media_espacial_df_corrientes["fecha"] = pd.to_datetime(
    cldamt_bajas_media_espacial_df_corrientes["fecha"])

# cldamt_medias
cldamt_medias_media_espacial_df_sudamerica = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_sudamerica.csv", index_col=0)
cldamt_medias_media_espacial_df_sudamerica["fecha"] = pd.to_datetime(
    cldamt_medias_media_espacial_df_sudamerica["fecha"])
cldamt_medias_media_espacial_df_region1 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_region1.csv", index_col=0)
cldamt_medias_media_espacial_df_region1["fecha"] = pd.to_datetime(
    cldamt_medias_media_espacial_df_region1["fecha"])
cldamt_medias_media_espacial_df_region2 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_region2.csv", index_col=0)
cldamt_medias_media_espacial_df_region2["fecha"] = pd.to_datetime(
    cldamt_medias_media_espacial_df_region2["fecha"])
cldamt_medias_media_espacial_df_corrientes = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_medias_media_espacial_df_corrientes.csv", index_col=0)
cldamt_medias_media_espacial_df_corrientes["fecha"] = pd.to_datetime(
    cldamt_medias_media_espacial_df_corrientes["fecha"])

# cldamt_altas
cldamt_altas_media_espacial_df_sudamerica = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_sudamerica.csv", index_col=0)
cldamt_altas_media_espacial_df_sudamerica["fecha"] = pd.to_datetime(
    cldamt_altas_media_espacial_df_sudamerica["fecha"])
cldamt_altas_media_espacial_df_region1 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_region1.csv", index_col=0)
cldamt_altas_media_espacial_df_region1["fecha"] = pd.to_datetime(
    cldamt_altas_media_espacial_df_region1["fecha"])
cldamt_altas_media_espacial_df_region2 = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_region2.csv", index_col=0)
cldamt_altas_media_espacial_df_region2["fecha"] = pd.to_datetime(
    cldamt_altas_media_espacial_df_region2["fecha"])
cldamt_altas_media_espacial_df_corrientes = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_altas_media_espacial_df_corrientes.csv", index_col=0)
cldamt_altas_media_espacial_df_corrientes["fecha"] = pd.to_datetime(
    cldamt_altas_media_espacial_df_corrientes["fecha"])


# %%

# separo por mes
def media_espacial_mes(media_espacial_df, numero_mes):
    media_espacial_df_mes = media_espacial_df[(
        pd.DatetimeIndex(media_espacial_df["fecha"]).month) == numero_mes]
    return(media_espacial_df_mes)
# %% corro para cldamt
# armo una lista con los df de media espacial mensual dde cldamt para cada region donde cada elemento tenga la serie de cada mes


# sudamerica
data_list_cldamt_media_espacial_mensuales_sudamerica = [None]*12
for i in range(0, 12):
    data_list_cldamt_media_espacial_mensuales_sudamerica[i] = media_espacial_mes(
        cldamt_media_espacial_df_sudamerica, i+1)

# region1
data_list_cldamt_media_espacial_mensuales_region1 = [None]*12
for i in range(0, 12):
    data_list_cldamt_media_espacial_mensuales_region1[i] = media_espacial_mes(
        cldamt_media_espacial_df_region1, i+1)

# region2
data_list_cldamt_media_espacial_mensuales_region2 = [None]*12
for i in range(0, 12):
    data_list_cldamt_media_espacial_mensuales_region2[i] = media_espacial_mes(
        cldamt_media_espacial_df_region2, i+1)

# corrientes
data_list_cldamt_media_espacial_mensuales_corrientes = [None]*12
for i in range(0, 12):
    data_list_cldamt_media_espacial_mensuales_corrientes[i] = media_espacial_mes(
        cldamt_media_espacial_df_corrientes, i+1)

# %% corro para cldamt_bajas
# armo una lista con los df de media espacial mensual dde cldamt para cada region donde cada elemento tenga la serie de cada mes

# sudamerica
data_list_cldamt_bajas_media_espacial_mensuales_sudamerica = [None]*12
for i in range(0, 12):
    data_list_cldamt_bajas_media_espacial_mensuales_sudamerica[i] = media_espacial_mes(
        cldamt_bajas_media_espacial_df_sudamerica, i+1)

# region1
data_list_cldamt_bajas_media_espacial_mensuales_region1 = [None]*12
for i in range(0, 12):
    data_list_cldamt_bajas_media_espacial_mensuales_region1[i] = media_espacial_mes(
        cldamt_bajas_media_espacial_df_region1, i+1)

# region2
data_list_cldamt_bajas_media_espacial_mensuales_region2 = [None]*12
for i in range(0, 12):
    data_list_cldamt_bajas_media_espacial_mensuales_region2[i] = media_espacial_mes(
        cldamt_bajas_media_espacial_df_region2, i+1)

# corrientes
data_list_cldamt_bajas_media_espacial_mensuales_corrientes = [None]*12
for i in range(0, 12):
    data_list_cldamt_bajas_media_espacial_mensuales_corrientes[i] = media_espacial_mes(
        cldamt_bajas_media_espacial_df_corrientes, i+1)

# %% corro para cldamt_medias
# armo una lista con los df de media espacial mensual dde cldamt para cada region donde cada elemento tenga la serie de cada mes

# sudamerica
data_list_cldamt_medias_media_espacial_mensuales_sudamerica = [None]*12
for i in range(0, 12):
    data_list_cldamt_medias_media_espacial_mensuales_sudamerica[i] = media_espacial_mes(
        cldamt_medias_media_espacial_df_sudamerica, i+1)

# region1
data_list_cldamt_medias_media_espacial_mensuales_region1 = [None]*12
for i in range(0, 12):
    data_list_cldamt_medias_media_espacial_mensuales_region1[i] = media_espacial_mes(
        cldamt_medias_media_espacial_df_region1, i+1)

# region2
data_list_cldamt_medias_media_espacial_mensuales_region2 = [None]*12
for i in range(0, 12):
    data_list_cldamt_medias_media_espacial_mensuales_region2[i] = media_espacial_mes(
        cldamt_medias_media_espacial_df_region2, i+1)

# corrientes
data_list_cldamt_medias_media_espacial_mensuales_corrientes = [None]*12
for i in range(0, 12):
    data_list_cldamt_medias_media_espacial_mensuales_corrientes[i] = media_espacial_mes(
        cldamt_medias_media_espacial_df_corrientes, i+1)

# %% corro para cldamt_altas
# armo una lista con los df de media espacial mensual dde cldamt para cada region donde cada elemento tenga la serie de cada mes

# sudamerica
data_list_cldamt_altas_media_espacial_mensuales_sudamerica = [None]*12
for i in range(0, 12):
    data_list_cldamt_altas_media_espacial_mensuales_sudamerica[i] = media_espacial_mes(
        cldamt_altas_media_espacial_df_sudamerica, i+1)

# region1
data_list_cldamt_altas_media_espacial_mensuales_region1 = [None]*12
for i in range(0, 12):
    data_list_cldamt_altas_media_espacial_mensuales_region1[i] = media_espacial_mes(
        cldamt_altas_media_espacial_df_region1, i+1)

# region2
data_list_cldamt_altas_media_espacial_mensuales_region2 = [None]*12
for i in range(0, 12):
    data_list_cldamt_altas_media_espacial_mensuales_region2[i] = media_espacial_mes(
        cldamt_altas_media_espacial_df_region2, i+1)

# corrientes
data_list_cldamt_altas_media_espacial_mensuales_corrientes = [None]*12
for i in range(0, 12):
    data_list_cldamt_altas_media_espacial_mensuales_corrientes[i] = media_espacial_mes(
        cldamt_altas_media_espacial_df_corrientes, i+1)

# %%
# separo por estacion


def media_espacial_estacion(media_espacial_df, numero_mes1, numero_mes2, numero_mes3, variable):
    media_espacial_df_meses_estacion = media_espacial_df[(((pd.DatetimeIndex(media_espacial_df["fecha"]).month) == numero_mes1) | (
        (pd.DatetimeIndex(media_espacial_df["fecha"]).month) == numero_mes2) | ((pd.DatetimeIndex(media_espacial_df["fecha"]).month) == numero_mes3))]
    media_espacial_df_meses_estacion.index = range(
        0, len(media_espacial_df_meses_estacion))
    media_espacial_df_media_estacion = pd.DataFrame(
        columns=["fecha", "Media_espacial_media_estacion"])
    for i in range(0, len(media_espacial_df_meses_estacion), 3):
        media_espacial_df_media_estacion = media_espacial_df_media_estacion.append({"fecha": pd.DatetimeIndex(media_espacial_df_meses_estacion["fecha"]).year[i], "Media_espacial_media_estacion": (
            (media_espacial_df_meses_estacion["Media_espacial_"+variable][i]+media_espacial_df_meses_estacion["Media_espacial_"+variable][i+1]+media_espacial_df_meses_estacion["Media_espacial_"+variable][i+2])/3)}, ignore_index=True)
    return(media_espacial_df_media_estacion)

# %% corro para cldamt
# armo una lista con los df de media espacial promedio estacional para cada anio dde cldamt para cada region donde cada elemento tenga la serie de cada estacion


# sudamerica
data_list_cldamt_media_espacial_estacional_sudamerica = [None]*4
data_list_cldamt_media_espacial_estacional_sudamerica[0] = media_espacial_estacion(
    cldamt_media_espacial_df_sudamerica, 12, 1, 2, "cldamt")
data_list_cldamt_media_espacial_estacional_sudamerica[1] = media_espacial_estacion(
    cldamt_media_espacial_df_sudamerica, 3, 4, 5, "cldamt")
data_list_cldamt_media_espacial_estacional_sudamerica[2] = media_espacial_estacion(
    cldamt_media_espacial_df_sudamerica, 6, 7, 8, "cldamt")
data_list_cldamt_media_espacial_estacional_sudamerica[3] = media_espacial_estacion(
    cldamt_media_espacial_df_sudamerica, 9, 10, 11, "cldamt")

# region1
data_list_cldamt_media_espacial_estacional_region1 = [None]*4
data_list_cldamt_media_espacial_estacional_region1[0] = media_espacial_estacion(
    cldamt_media_espacial_df_region1, 12, 1, 2, "cldamt")
data_list_cldamt_media_espacial_estacional_region1[1] = media_espacial_estacion(
    cldamt_media_espacial_df_region1, 3, 4, 5, "cldamt")
data_list_cldamt_media_espacial_estacional_region1[2] = media_espacial_estacion(
    cldamt_media_espacial_df_region1, 6, 7, 8, "cldamt")
data_list_cldamt_media_espacial_estacional_region1[3] = media_espacial_estacion(
    cldamt_media_espacial_df_region1, 9, 10, 11, "cldamt")

# region2
data_list_cldamt_media_espacial_estacional_region2 = [None]*4
data_list_cldamt_media_espacial_estacional_region2[0] = media_espacial_estacion(
    cldamt_media_espacial_df_region2, 12, 1, 2, "cldamt")
data_list_cldamt_media_espacial_estacional_region2[1] = media_espacial_estacion(
    cldamt_media_espacial_df_region2, 3, 4, 5, "cldamt")
data_list_cldamt_media_espacial_estacional_region2[2] = media_espacial_estacion(
    cldamt_media_espacial_df_region2, 6, 7, 8, "cldamt")
data_list_cldamt_media_espacial_estacional_region2[3] = media_espacial_estacion(
    cldamt_media_espacial_df_region2, 9, 10, 11, "cldamt")

# corrientes
data_list_cldamt_media_espacial_estacional_corrientes = [None]*4
data_list_cldamt_media_espacial_estacional_corrientes[0] = media_espacial_estacion(
    cldamt_media_espacial_df_corrientes, 12, 1, 2, "cldamt")
data_list_cldamt_media_espacial_estacional_corrientes[1] = media_espacial_estacion(
    cldamt_media_espacial_df_corrientes, 3, 4, 5, "cldamt")
data_list_cldamt_media_espacial_estacional_corrientes[2] = media_espacial_estacion(
    cldamt_media_espacial_df_corrientes, 6, 7, 8, "cldamt")
data_list_cldamt_media_espacial_estacional_corrientes[3] = media_espacial_estacion(
    cldamt_media_espacial_df_corrientes, 9, 10, 11, "cldamt")

# %% corro para cldamt bajas
# armo una lista con los df de media espacial promedio estacional para cada anio de cldamt bajas para cada region donde cada elemento tenga la serie de cada estacion

# sudamerica
data_list_cldamt_bajas_media_espacial_estacional_sudamerica = [None]*4
data_list_cldamt_bajas_media_espacial_estacional_sudamerica[0] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_sudamerica, 12, 1, 2, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_sudamerica[1] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_sudamerica, 3, 4, 5, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_sudamerica[2] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_sudamerica, 6, 7, 8, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_sudamerica[3] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_sudamerica, 9, 10, 11, "cldamt_bajas")

# region1
data_list_cldamt_bajas_media_espacial_estacional_region1 = [None]*4
data_list_cldamt_bajas_media_espacial_estacional_region1[0] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region1, 12, 1, 2, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region1[1] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region1, 3, 4, 5, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region1[2] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region1, 6, 7, 8, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region1[3] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region1, 9, 10, 11, "cldamt_bajas")

# region2
data_list_cldamt_bajas_media_espacial_estacional_region2 = [None]*4
data_list_cldamt_bajas_media_espacial_estacional_region2[0] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region2, 12, 1, 2, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region2[1] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region2, 3, 4, 5, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region2[2] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region2, 6, 7, 8, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_region2[3] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_region2, 9, 10, 11, "cldamt_bajas")

# corrientes
data_list_cldamt_bajas_media_espacial_estacional_corrientes = [None]*4
data_list_cldamt_bajas_media_espacial_estacional_corrientes[0] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_corrientes, 12, 1, 2, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_corrientes[1] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_corrientes, 3, 4, 5, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_corrientes[2] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_corrientes, 6, 7, 8, "cldamt_bajas")
data_list_cldamt_bajas_media_espacial_estacional_corrientes[3] = media_espacial_estacion(
    cldamt_bajas_media_espacial_df_corrientes, 9, 10, 11, "cldamt_bajas")

# %% corro para cldamt medias
# armo una lista con los df de media espacial promedio estacional para cada anio de cldamt bajas para cada region donde cada elemento tenga la serie de cada estacion

# sudamerica
data_list_cldamt_medias_media_espacial_estacional_sudamerica = [None]*4
data_list_cldamt_medias_media_espacial_estacional_sudamerica[0] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_sudamerica, 12, 1, 2, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_sudamerica[1] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_sudamerica, 3, 4, 5, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_sudamerica[2] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_sudamerica, 6, 7, 8, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_sudamerica[3] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_sudamerica, 9, 10, 11, "cldamt_medias")

# region1
data_list_cldamt_medias_media_espacial_estacional_region1 = [None]*4
data_list_cldamt_medias_media_espacial_estacional_region1[0] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region1, 12, 1, 2, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region1[1] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region1, 3, 4, 5, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region1[2] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region1, 6, 7, 8, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region1[3] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region1, 9, 10, 11, "cldamt_medias")

# region2
data_list_cldamt_medias_media_espacial_estacional_region2 = [None]*4
data_list_cldamt_medias_media_espacial_estacional_region2[0] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region2, 12, 1, 2, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region2[1] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region2, 3, 4, 5, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region2[2] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region2, 6, 7, 8, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_region2[3] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_region2, 9, 10, 11, "cldamt_medias")

# corrientes
data_list_cldamt_medias_media_espacial_estacional_corrientes = [None]*4
data_list_cldamt_medias_media_espacial_estacional_corrientes[0] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_corrientes, 12, 1, 2, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_corrientes[1] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_corrientes, 3, 4, 5, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_corrientes[2] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_corrientes, 6, 7, 8, "cldamt_medias")
data_list_cldamt_medias_media_espacial_estacional_corrientes[3] = media_espacial_estacion(
    cldamt_medias_media_espacial_df_corrientes, 9, 10, 11, "cldamt_medias")

# %% corro para cldamt altas
# armo una lista con los df de media espacial promedio estacional para cada anio de cldamt altas para cada region donde cada elemento tenga la serie de cada estacion

# sudamerica
data_list_cldamt_altas_media_espacial_estacional_sudamerica = [None]*4
data_list_cldamt_altas_media_espacial_estacional_sudamerica[0] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_sudamerica, 12, 1, 2, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_sudamerica[1] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_sudamerica, 3, 4, 5, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_sudamerica[2] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_sudamerica, 6, 7, 8, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_sudamerica[3] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_sudamerica, 9, 10, 11, "cldamt_altas")

# region1
data_list_cldamt_altas_media_espacial_estacional_region1 = [None]*4
data_list_cldamt_altas_media_espacial_estacional_region1[0] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region1, 12, 1, 2, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region1[1] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region1, 3, 4, 5, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region1[2] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region1, 6, 7, 8, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region1[3] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region1, 9, 10, 11, "cldamt_altas")

# region2
data_list_cldamt_altas_media_espacial_estacional_region2 = [None]*4
data_list_cldamt_altas_media_espacial_estacional_region2[0] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region2, 12, 1, 2, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region2[1] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region2, 3, 4, 5, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region2[2] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region2, 6, 7, 8, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_region2[3] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_region2, 9, 10, 11, "cldamt_altas")

# corrientes
data_list_cldamt_altas_media_espacial_estacional_corrientes = [None]*4
data_list_cldamt_altas_media_espacial_estacional_corrientes[0] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_corrientes, 12, 1, 2, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_corrientes[1] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_corrientes, 3, 4, 5, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_corrientes[2] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_corrientes, 6, 7, 8, "cldamt_altas")
data_list_cldamt_altas_media_espacial_estacional_corrientes[3] = media_espacial_estacion(
    cldamt_altas_media_espacial_df_corrientes, 9, 10, 11, "cldamt_altas")


# %%
"""
grafico las series temporales
"""
# armo funcion para graficar serie completa
# armo funcion para graficar series de cada estacion (en un mismo plot)
# armo funcion para graficar las series por cada anio (en un mismo plot)

# anual


def serie_periodo_completo(data_frame_entrada, variable, region, ruta_salida):
    fig1, ax = plt.subplots(figsize=[10, 5], dpi=200)
    plt.plot(data_frame_entrada["fecha"],
             data_frame_entrada["Media_espacial_"+variable], color="indigo")
    ax.tick_params(axis='x', direction='out', bottom=True,
                   labelrotation=45, labelsize=10, pad=1.5)
    # ax.set_ylim(20,80)
    ax.set_xlabel("Fecha", size=10)
    ax.set_ylabel(variable+" %", size=10)
    ax.grid()
    plt.title(variable+" Media mensual media "+region + " (serie completa)")
    nombre = variable+"_"+"media_espacial_mensual_" + \
        region+"_"+"(serie completa)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %%corro para cldamt
# corro la funcion para graficar la serie de cada region
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"
serie_periodo_completo(cldamt_media_espacial_df_sudamerica,
                       "cldamt", "Sudamérica", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_region1,
                       "cldamt", "Región 1", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_region2,
                       "cldamt", "Región 2", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_corrientes,
                       "cldamt", "Corrientes", ruta_salida)

# %%corro para cldamt_bajas
# corro la funcion para graficar la serie de cada region
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_series/"
serie_periodo_completo(cldamt_bajas_media_espacial_df_sudamerica,
                       "cldamt_bajas", "Sudamérica", ruta_salida)
serie_periodo_completo(cldamt_bajas_media_espacial_df_region1,
                       "cldamt_bajas", "Región 1", ruta_salida)
serie_periodo_completo(cldamt_bajas_media_espacial_df_region2,
                       "cldamt_bajas", "Región 2", ruta_salida)
serie_periodo_completo(cldamt_bajas_media_espacial_df_corrientes,
                       "cldamt_bajas", "Corrientes", ruta_salida)

# %%corro para cldamt_medias
# corro la funcion para graficar la serie de cada region
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_series/"
serie_periodo_completo(cldamt_medias_media_espacial_df_sudamerica,
                       "cldamt_medias", "Sudamérica", ruta_salida)
serie_periodo_completo(cldamt_medias_media_espacial_df_region1,
                       "cldamt_medias", "Región 1", ruta_salida)
serie_periodo_completo(cldamt_medias_media_espacial_df_region2,
                       "cldamt_medias", "Región 2", ruta_salida)
serie_periodo_completo(cldamt_medias_media_espacial_df_corrientes,
                       "cldamt_medias", "Corrientes", ruta_salida)

# %%corro para cldamt_altas
# corro la funcion para graficar la serie de cada region
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_series/"
serie_periodo_completo(cldamt_altas_media_espacial_df_sudamerica,
                       "cldamt_altas", "Sudamérica", ruta_salida)
serie_periodo_completo(cldamt_altas_media_espacial_df_region1,
                       "cldamt_altas", "Región 1", ruta_salida)
serie_periodo_completo(cldamt_altas_media_espacial_df_region2,
                       "cldamt_altas", "Región 2", ruta_salida)
serie_periodo_completo(cldamt_altas_media_espacial_df_corrientes,
                       "cldamt_altas", "Corrientes", ruta_salida)


# %%
# mensual

def serie_mensual(lista, variable, region, ymin, ymax, ruta_salida):

    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1, ax = plt.subplots(4, 3, figsize=[12, 10], dpi=200)
    fig1.suptitle(variable+" Media mensual media " +
                  region + " (meses)", size=18)

    meses = [["Enero", "Febrero", "Marzo"], ["Abril", "Mayo", "Junio"], [
        "Julio", "Agosto", "Septiembre"], ["Octubre", "Noviembre", "Diciembre"]]
    for j in range(0, 4):
        for i in range(0, 3):
            ax[j, i].plot(lista[i+3*j]["fecha"], lista[i+3*j]
                          ["Media_espacial_"+variable], color="indigo")
            ax[j, i].tick_params(
                axis='x', direction='out', bottom=True, labelrotation=45, labelsize=10, pad=1.5)
            ax[j, i].set_ylim(ymin, ymax)
            ax[j, i].set_xlabel("Fecha", size=10)
            ax[j, i].set_ylabel(variable+" %", size=10)
            ax[j, i].grid()
            ax[j, i].set_title(meses[j][i])

    fig1.tight_layout()
    nombre = variable+"_"+"media_espacial_mensual_"+region+"_"+"(meses)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %% corro para cldamt
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"
serie_mensual(data_list_cldamt_media_espacial_mensuales_sudamerica,
              "cldamt", "Sudamérica", 60, 80, ruta_salida)
serie_mensual(data_list_cldamt_media_espacial_mensuales_region1,
              "cldamt", "Región 1", 50, 80, ruta_salida)
serie_mensual(data_list_cldamt_media_espacial_mensuales_region2,
              "cldamt", "Región 2", 30, 80, ruta_salida)
serie_mensual(data_list_cldamt_media_espacial_mensuales_corrientes,
              "cldamt", "Corrientes", 20, 80, ruta_salida)
# %% corro para cldamt_bajas
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_series/"
serie_mensual(data_list_cldamt_bajas_media_espacial_mensuales_sudamerica,
              "cldamt_bajas", "Sudamérica", 0, 40, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_bajas_media_espacial_mensuales_region1,
              "cldamt_bajas", "Región 1", 0, 40, ruta_salida)
serie_mensual(data_list_cldamt_bajas_media_espacial_mensuales_region2,
              "cldamt_bajas", "Región 2", 0, 40, ruta_salida)
serie_mensual(data_list_cldamt_bajas_media_espacial_mensuales_corrientes,
              "cldamt_bajas", "Corrientes", 0, 40, ruta_salida)
# %% corro para cldamt_medias
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_series/"
serie_mensual(data_list_cldamt_medias_media_espacial_mensuales_sudamerica,
              "cldamt_medias", "Sudamérica", 0, 40, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_medias_media_espacial_mensuales_region1,
              "cldamt_medias", "Región 1", 0, 40, ruta_salida)
serie_mensual(data_list_cldamt_medias_media_espacial_mensuales_region2,
              "cldamt_medias", "Región 2", 0, 40, ruta_salida)
serie_mensual(data_list_cldamt_medias_media_espacial_mensuales_corrientes,
              "cldamt_medias", "Corrientes", 0, 40, ruta_salida)

# %% corro para cldamt_altas
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_series/"
serie_mensual(data_list_cldamt_altas_media_espacial_mensuales_sudamerica,
              "cldamt_altas", "Sudamérica", 10, 50, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_altas_media_espacial_mensuales_region1,
              "cldamt_altas", "Región 1", 10, 50, ruta_salida)
serie_mensual(data_list_cldamt_altas_media_espacial_mensuales_region2,
              "cldamt_altas", "Región 2", 10, 50, ruta_salida)
serie_mensual(data_list_cldamt_altas_media_espacial_mensuales_corrientes,
              "cldamt_altas", "Corrientes", 10, 50, ruta_salida)


# %%
# trimestral


def serie_trimestral(lista, variable, region, ymin, ymax, ruta_salida):

    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1, ax = plt.subplots(2, 2, figsize=[12, 10], dpi=200)
    fig1.suptitle(variable+" Media trimestral media " +
                  region + " (estaciones)", size=18)

    estacion = [["DEF", "MAM"], ["JJA", "SON"]]
    for j in range(0, 2):
        for i in range(0, 2):
            ax[j, i].plot(lista[i+2*j]["fecha"], lista[i+2*j]
                          ["Media_espacial_media_estacion"], color="indigo")
            ax[j, i].tick_params(
                axis='x', direction='out', bottom=True, labelrotation=45, labelsize=10, pad=1.5)
            ax[j, i].set_ylim(ymin, ymax)
            ax[j, i].set_xlabel("Fecha", size=10)
            ax[j, i].set_ylabel(variable+" %", size=10)
            ax[j, i].grid()
            ax[j, i].set_title(estacion[j][i])

    fig1.tight_layout()
    nombre = variable+"_"+"media_espacial_trimestral_" + \
        region+"_"+"(estaciones)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %% corro para cldamt
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_series/"
serie_trimestral(data_list_cldamt_media_espacial_estacional_sudamerica,
                 "cldamt", "Sudamérica", 60, 80, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region1,
                 "cldamt", "Región 1", 50, 80, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region2,
                 "cldamt", "Región 2", 40, 80, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_corrientes,
                 "cldamt", "Corrientes", 40, 80, ruta_salida)

# %% corro para cldamt_bajas
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_series/"
serie_trimestral(data_list_cldamt_bajas_media_espacial_estacional_sudamerica,
                 "cldamt_bajas", "Sudamérica", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_bajas_media_espacial_estacional_region1,
                 "cldamt_bajas", "Región 1", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_bajas_media_espacial_estacional_region2,
                 "cldamt_bajas", "Región 2", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_bajas_media_espacial_estacional_corrientes,
                 "cldamt_bajas", "Corrientes", 0, 40, ruta_salida)

# %% corro para cldamt_medias
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_medias_series/"
serie_trimestral(data_list_cldamt_medias_media_espacial_estacional_sudamerica,
                 "cldamt_medias", "Sudamérica", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_medias_media_espacial_estacional_region1,
                 "cldamt_medias", "Región 1", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_medias_media_espacial_estacional_region2,
                 "cldamt_medias", "Región 2", 0, 40, ruta_salida)
serie_trimestral(data_list_cldamt_medias_media_espacial_estacional_corrientes,
                 "cldamt_medias", "Corrientes", 0, 40, ruta_salida)

# %% corro para cldamt_altas
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_altas_series/"
serie_trimestral(data_list_cldamt_altas_media_espacial_estacional_sudamerica,
                 "cldamt_altas", "Sudamérica", 10, 50, ruta_salida)
serie_trimestral(data_list_cldamt_altas_media_espacial_estacional_region1,
                 "cldamt_altas", "Región 1", 10, 50, ruta_salida)
serie_trimestral(data_list_cldamt_altas_media_espacial_estacional_region2,
                 "cldamt_altas", "Región 2", 10, 50, ruta_salida)
serie_trimestral(data_list_cldamt_altas_media_espacial_estacional_corrientes,
                 "cldamt_altas", "Corrientes", 10, 50, ruta_salida)


# %%
"""
Veo tendencia y test mann kendall
"""

# cldamt
cldamt_media_espacial_df_sudamerica = pd.read_csv(
    "Documentos/Doctorado/datos/nubosidad/cldamt_media_espacial_df_sudamerica.csv", index_col=0)

data = cldamt_media_espacial_df_sudamerica["Media_espacial_cldamt"]
mk.original_test(data, alpha=0.05)

# https://stackoverflow.com/questions/6148207/linear-regression-with-matplotlib-numpy

data = np.array(data)
fecha = cldamt_media_espacial_df_sudamerica["fecha"]
fecha_list = np.arange(0, len(fecha), 1)
coef = np.polyfit(fecha_list, data, 1)
poly1d_fn = np.poly1d(coef)
poly1d_fn(fecha_list)

plt.plot(fecha_list, data, '-r', fecha_list, poly1d_fn(fecha_list), '--k')

media = coef[1]
# esto es mensual. Para hacerlo decadal multiplico por 12 y por 10
tendencia = coef[0]
tendencia_decadal = coef[0]*12*10

desvio = np.std(data)

# %% lo hago con tendencia y testeandola con mann kendall

"""
grafico las series temporales completas para todos los tipos de nube juntos 
"""
# anual


def serie_periodo_completo(data_frame_entrada_cldamt, data_frame_entrada_bajas, data_frame_entrada_medias, data_frame_entrada_altas, variable, region, ruta_salida):
    fig1, ax = plt.subplots(figsize=[10, 6], dpi=200)
    plt.plot(data_frame_entrada_cldamt["fecha"],
             data_frame_entrada_cldamt["Media_espacial_"+variable], color="k", alpha=0.7)
    plt.plot(data_frame_entrada_bajas["fecha"],
             data_frame_entrada_bajas["Media_espacial_"+variable+"_bajas"], color="teal", alpha=0.7)
    plt.plot(data_frame_entrada_medias["fecha"],
             data_frame_entrada_medias["Media_espacial_"+variable+"_medias"], color="purple", alpha=0.7)
    plt.plot(data_frame_entrada_altas["fecha"], data_frame_entrada_altas["Media_espacial_" +
             variable+"_altas"], color="crimson", alpha=0.7)

    # agrego lineas de tendencia, si es significativa con un 95% (test mann kendall) lo hago con linea llena y si no con linea intermitente
    coef1 = np.polyfit(np.arange(0, len(data_frame_entrada_cldamt["fecha"]), 1), np.array(
        data_frame_entrada_cldamt["Media_espacial_"+variable]), 1)
    linear_fit_1 = np.poly1d(coef1)
    if (abs(mk.original_test(data_frame_entrada_cldamt["Media_espacial_"+variable], alpha=0.05)[3]) > 1.96):
        plt.plot(data_frame_entrada_cldamt["fecha"], linear_fit_1(
            np.arange(0, len(data_frame_entrada_cldamt["fecha"]), 1)), color="k", ls="-")
    else:
        plt.plot(data_frame_entrada_cldamt["fecha"], linear_fit_1(
            np.arange(0, len(data_frame_entrada_cldamt["fecha"]), 1)), color="k", ls=":")

    coef2 = np.polyfit(np.arange(0, len(data_frame_entrada_bajas["fecha"]), 1), np.array(
        data_frame_entrada_bajas["Media_espacial_"+variable+"_bajas"]), 1)
    linear_fit_2 = np.poly1d(coef2)
    if (abs(mk.original_test(data_frame_entrada_bajas["Media_espacial_"+variable+"_bajas"], alpha=0.05)[3]) > 1.96):
        plt.plot(data_frame_entrada_bajas["fecha"], linear_fit_2(
            np.arange(0, len(data_frame_entrada_bajas["fecha"]), 1)), color="teal", ls="-")
    else:
        plt.plot(data_frame_entrada_bajas["fecha"], linear_fit_2(
            np.arange(0, len(data_frame_entrada_bajas["fecha"]), 1)), color="teal", ls=":")

    coef3 = np.polyfit(np.arange(0, len(data_frame_entrada_medias["fecha"]), 1), np.array(
        data_frame_entrada_medias["Media_espacial_"+variable+"_medias"]), 1)
    linear_fit_3 = np.poly1d(coef3)
    if (abs(mk.original_test(data_frame_entrada_medias["Media_espacial_"+variable+"_medias"], alpha=0.05)[3]) > 1.96):
        plt.plot(data_frame_entrada_medias["fecha"], linear_fit_3(np.arange(
            0, len(data_frame_entrada_medias["fecha"]), 1)), color="purple", ls="-")
    else:
        plt.plot(data_frame_entrada_medias["fecha"], linear_fit_3(np.arange(
            0, len(data_frame_entrada_medias["fecha"]), 1)), color="purple", ls=":")

    coef4 = np.polyfit(np.arange(0, len(data_frame_entrada_altas["fecha"]), 1), np.array(
        data_frame_entrada_altas["Media_espacial_"+variable+"_altas"]), 1)
    linear_fit_4 = np.poly1d(coef4)
    if (abs(mk.original_test(data_frame_entrada_altas["Media_espacial_"+variable+"_altas"], alpha=0.05)[3]) > 1.96):
        plt.plot(data_frame_entrada_altas["fecha"], linear_fit_4(np.arange(
            0, len(data_frame_entrada_altas["fecha"]), 1)), color="crimson", ls="-")
    else:
        plt.plot(data_frame_entrada_altas["fecha"], linear_fit_4(np.arange(
            0, len(data_frame_entrada_altas["fecha"]), 1)), color="crimson", ls="--")

    media1 = np.round(np.nanmean(
        data_frame_entrada_cldamt["Media_espacial_"+variable]), 1)
    desvio1 = np.round(
        np.nanstd(data_frame_entrada_cldamt["Media_espacial_"+variable]), 2)
    tendencia1 = np.round(coef1[0]*12*10, 2)  # decadal

    media2 = np.round(np.nanmean(
        data_frame_entrada_bajas["Media_espacial_"+variable+"_bajas"]), 1)
    desvio2 = np.round(
        np.nanstd(data_frame_entrada_bajas["Media_espacial_"+variable+"_bajas"]), 2)
    tendencia2 = np.round(coef2[0]*12*10, 2)  # decadal

    media3 = np.round(np.nanmean(
        data_frame_entrada_medias["Media_espacial_"+variable+"_medias"]), 1)
    desvio3 = np.round(
        np.nanstd(data_frame_entrada_medias["Media_espacial_"+variable+"_medias"]), 2)
    tendencia3 = np.round(coef3[0]*12*10, 2)  # decadal

    media4 = np.round(np.nanmean(
        data_frame_entrada_altas["Media_espacial_"+variable+"_altas"]), 1)
    desvio4 = np.round(
        np.nanstd(data_frame_entrada_altas["Media_espacial_"+variable+"_altas"]), 2)
    tendencia4 = np.round(coef4[0]*12*10, 2)  # decadal

    plt.text(5200, 85, "Media (%) \nDesvío (%) \nTendencia decadal (%/dec)",
             color="k", ha="left", backgroundcolor="white")
    plt.text(9500, 85, str(media1)+" \n"+str(desvio1)+" \n" +
             str(tendencia1), color="k", ha="left", backgroundcolor="white")
    plt.text(10500, 85, str(media2)+" \n"+str(desvio2)+" \n" +
             str(tendencia2), color="teal", ha="left", backgroundcolor="white")
    plt.text(11500, 85, str(media3)+" \n"+str(desvio3)+" \n" +
             str(tendencia3), color="purple", ha="left", backgroundcolor="white")
    plt.text(12500, 85, str(media4)+" \n"+str(desvio4)+" \n" +
             str(tendencia4), color="crimson", ha="left", backgroundcolor="white")

    ax.tick_params(axis='x', direction='out', bottom=True,
                   labelrotation=25, labelsize=10, pad=1.5)
    ax.set_ylim(0, 100)
    # plt.xticks(data_frame_entrada_cldamt["fecha"][::24])
    ax.set_xlabel("Fecha", size=10)
    ax.set_ylabel(variable+" %", size=10)
    ax.grid()
    plt.title(variable+" Media mensual media "+region + " (serie completa)")
    plt.legend(["cldamt", "cldamt bajas", "cldamt medias", "cldamt altas"])
    nombre = variable+"_"+"media_espacial_mensual_" + \
        region+"_"+"(serie completa)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %%
# lo corro
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_medias_altas_series/"
serie_periodo_completo(cldamt_media_espacial_df_sudamerica, cldamt_bajas_media_espacial_df_sudamerica,
                       cldamt_medias_media_espacial_df_sudamerica, cldamt_altas_media_espacial_df_sudamerica, "cldamt", "Sudamérica", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_region1, cldamt_bajas_media_espacial_df_region1,
                       cldamt_medias_media_espacial_df_region1, cldamt_altas_media_espacial_df_region1, "cldamt", "Región 1", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_region2, cldamt_bajas_media_espacial_df_region2,
                       cldamt_medias_media_espacial_df_region2, cldamt_altas_media_espacial_df_region2, "cldamt", "Región 2", ruta_salida)
serie_periodo_completo(cldamt_media_espacial_df_corrientes, cldamt_bajas_media_espacial_df_corrientes,
                       cldamt_medias_media_espacial_df_corrientes, cldamt_altas_media_espacial_df_corrientes, "cldamt", "Corrientes", ruta_salida)

# %%
"""
Lo hago para cldamt, bajas, medias y altas junto  VER DE PONER SOLO LOS NUMEROS DE LOS ANIOS EN LOS EJES 
"""
# mensual


def serie_mensual(lista_cldamt, lista_bajas, lista_medias, lista_altas, variable, region, ymin, ymax, ruta_salida):

    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1, ax = plt.subplots(4, 3, figsize=[12, 10], dpi=200)
    fig1.suptitle(variable+" Media mensual media " +
                  region + " (meses)", size=18)

    meses = [["Enero", "Febrero", "Marzo"], ["Abril", "Mayo", "Junio"], [
        "Julio", "Agosto", "Septiembre"], ["Octubre", "Noviembre", "Diciembre"]]
    for j in range(0, 4):
        for i in range(0, 3):
            ax[j, i].plot(lista_cldamt[i+3*j]["fecha"], lista_cldamt[i+3*j]
                          ["Media_espacial_"+variable], color="k", alpha=0.7)
            ax[j, i].plot(lista_bajas[i+3*j]["fecha"], lista_bajas[i+3*j]
                          ["Media_espacial_"+variable+"_bajas"], color="teal", alpha=0.7)
            ax[j, i].plot(lista_medias[i+3*j]["fecha"], lista_medias[i+3*j]
                          ["Media_espacial_"+variable+"_medias"], color="purple", alpha=0.7)
            ax[j, i].plot(lista_altas[i+3*j]["fecha"], lista_altas[i+3*j]
                          ["Media_espacial_"+variable+"_altas"], color="crimson", alpha=0.7)

            # agrego lineas de tendencia, si es significativa con un 95% (test mann kendall) lo hago con linea llena y si no con linea intermitente
            coef1 = np.polyfit(np.arange(0, len(lista_cldamt[i+3*j]["fecha"]), 1), np.array(
                lista_cldamt[i+3*j]["Media_espacial_"+variable]), 1)
            linear_fit_1 = np.poly1d(coef1)
            if (abs(mk.original_test(lista_cldamt[i+3*j]["Media_espacial_"+variable], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_cldamt[i+3*j]["fecha"], linear_fit_1(
                    np.arange(0, len(lista_cldamt[i+3*j]["fecha"]), 1)), color="k", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_cldamt[i+3*j]["fecha"], linear_fit_1(
                    np.arange(0, len(lista_cldamt[i+3*j]["fecha"]), 1)), color="k", ls=":", lw=0.9)

            coef2 = np.polyfit(np.arange(0, len(lista_bajas[i+3*j]["fecha"]), 1), np.array(
                lista_bajas[i+3*j]["Media_espacial_"+variable+"_bajas"]), 1)
            linear_fit_2 = np.poly1d(coef2)
            if (abs(mk.original_test(lista_bajas[i+3*j]["Media_espacial_"+variable+"_bajas"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_bajas[i+3*j]["fecha"], linear_fit_2(np.arange(
                    0, len(lista_bajas[i+3*j]["fecha"]), 1)), color="teal", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_bajas[i+3*j]["fecha"], linear_fit_2(np.arange(
                    0, len(lista_bajas[i+3*j]["fecha"]), 1)), color="teal", ls=":", lw=0.9)

            coef3 = np.polyfit(np.arange(0, len(lista_medias[i+3*j]["fecha"]), 1), np.array(
                lista_medias[i+3*j]["Media_espacial_"+variable+"_medias"]), 1)
            linear_fit_3 = np.poly1d(coef3)
            if (abs(mk.original_test(lista_medias[i+3*j]["Media_espacial_"+variable+"_medias"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_medias[i+3*j]["fecha"], linear_fit_3(np.arange(
                    0, len(lista_medias[i+3*j]["fecha"]), 1)), color="purple", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_medias[i+3*j]["fecha"], linear_fit_3(np.arange(
                    0, len(lista_medias[i+3*j]["fecha"]), 1)), color="purple", ls=":", lw=0.9)

            coef4 = np.polyfit(np.arange(0, len(lista_altas[i+3*j]["fecha"]), 1), np.array(
                lista_altas[i+3*j]["Media_espacial_"+variable+"_altas"]), 1)
            linear_fit_4 = np.poly1d(coef4)
            if (abs(mk.original_test(lista_altas[i+3*j]["Media_espacial_"+variable+"_altas"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_altas[i+3*j]["fecha"], linear_fit_4(np.arange(
                    0, len(lista_altas[i+3*j]["fecha"]), 1)), color="crimson", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_altas[i+3*j]["fecha"], linear_fit_4(np.arange(
                    0, len(lista_altas[i+3*j]["fecha"]), 1)), color="crimson", ls=":", lw=0.9)

            ax[j, i].tick_params(
                axis='x', direction='out', bottom=True, labelrotation=45, labelsize=10, pad=1.5)
            ax[j, i].set_ylim(ymin, ymax)
            ax[j, i].set_xlabel("Fecha", size=10)
            ax[j, i].set_ylabel(variable+" %", size=10)
            ax[j, i].grid()
            ax[j, i].set_title(meses[j][i])
            ax[j, i].text(8300, 75, str(np.round(np.nanmean(lista_cldamt[i+3*j]["Media_espacial_"+variable]), 1))+" \n"+str(np.round(np.nanstd(
                lista_cldamt[i+3*j]["Media_espacial_"+variable]), 2))+" \n"+np.str(np.round(coef1[0]*10, 2)), fontsize="x-small", c="k")  # , backgroundcolor="white")
            ax[j, i].text(9800, 75, str(np.round(np.nanmean(lista_bajas[i+3*j]["Media_espacial_"+variable+"_bajas"]), 1))+" \n"+str(np.round(np.nanstd(lista_bajas[i+3*j]
                          ["Media_espacial_"+variable+"_bajas"]), 2))+" \n"+np.str(np.round(coef2[0]*10, 2)), fontsize="x-small", c="teal")  # , backgroundcolor="white")
            ax[j, i].text(11300, 75, str(np.round(np.nanmean(lista_medias[i+3*j]["Media_espacial_"+variable+"_medias"]), 1))+" \n"+str(np.round(np.nanstd(lista_medias[i+3*j]
                          ["Media_espacial_"+variable+"_medias"]), 2))+" \n"+np.str(np.round(coef3[0]*10, 2)), fontsize="x-small", c="purple")  # , backgroundcolor="white")
            ax[j, i].text(12800, 75, str(np.round(np.nanmean(lista_altas[i+3*j]["Media_espacial_"+variable+"_altas"]), 1))+" \n"+str(np.round(np.nanstd(lista_altas[i+3*j]
                          ["Media_espacial_"+variable+"_altas"]), 2))+" \n"+np.str(np.round(coef4[0]*10, 2)), fontsize="x-small", c="crimson")  # , backgroundcolor="white")

    # plt.xticks(lista_cldamt[0]["fecha"][::1])
    fig1.legend(["cldamt", "cldamt bajas", "cldamt medias", "cldamt altas"],
                loc='lower left', ncol=4, bbox_to_anchor=(0.25, 0.00))
    # fig1.tight_layout
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.3,
                        hspace=0.7)
    nombre = variable+"_"+"media_espacial_mensual_"+region+"_"+"(meses)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %%
"""
Corro funcion
"""
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_medias_altas_series/"
serie_mensual(data_list_cldamt_media_espacial_mensuales_sudamerica, data_list_cldamt_bajas_media_espacial_mensuales_sudamerica,
              data_list_cldamt_medias_media_espacial_mensuales_sudamerica, data_list_cldamt_altas_media_espacial_mensuales_sudamerica, "cldamt", "Sudamérica", 0, 100, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_media_espacial_mensuales_region1, data_list_cldamt_bajas_media_espacial_mensuales_region1,
              data_list_cldamt_medias_media_espacial_mensuales_region1, data_list_cldamt_altas_media_espacial_mensuales_region1, "cldamt", "Región 1", 0, 100, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_media_espacial_mensuales_region2, data_list_cldamt_bajas_media_espacial_mensuales_region2,
              data_list_cldamt_medias_media_espacial_mensuales_region2, data_list_cldamt_altas_media_espacial_mensuales_region2, "cldamt", "Región 2", 0, 100, ruta_salida)  # 60,80
serie_mensual(data_list_cldamt_media_espacial_mensuales_corrientes, data_list_cldamt_bajas_media_espacial_mensuales_corrientes,
              data_list_cldamt_medias_media_espacial_mensuales_corrientes, data_list_cldamt_altas_media_espacial_mensuales_corrientes, "cldamt", "Corrientes", 0, 100, ruta_salida)  # 60,80

# %% armo funcion con cldamt, bajas, medias y altas juntas

# trimestral


def serie_trimestral(lista_cldamt, lista_bajas, lista_medias, lista_altas, variable, region, ymin, ymax, ruta_salida):

    # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/subplots_demo.html
    fig1, ax = plt.subplots(2, 2, figsize=[11, 11], dpi=200)
    fig1.suptitle(variable+" Media trimestral media " +
                  region + " (estaciones)", size=18)

    estacion = [["DEF", "MAM"], ["JJA", "SON"]]
    for j in range(0, 2):
        for i in range(0, 2):
            ax[j, i].plot(lista_cldamt[i+2*j]["fecha"], lista_cldamt[i+2*j]
                          ["Media_espacial_media_estacion"], color="k", alpha=0.7)
            ax[j, i].plot(lista_bajas[i+2*j]["fecha"], lista_bajas[i+2*j]
                          ["Media_espacial_media_estacion"], color="teal", alpha=0.7)
            ax[j, i].plot(lista_medias[i+2*j]["fecha"], lista_medias[i+2*j]
                          ["Media_espacial_media_estacion"], color="purple", alpha=0.7)
            ax[j, i].plot(lista_altas[i+2*j]["fecha"], lista_altas[i+2*j]
                          ["Media_espacial_media_estacion"], color="crimson", alpha=0.7)

            # agrego lineas de tendencia, si es significativa con un 95% (test mann kendall) lo hago con linea llena y si no con linea intermitente
            coef1 = np.polyfit(np.arange(0, len(lista_cldamt[i+2*j]["fecha"]), 1), np.array(
                lista_cldamt[i+2*j]["Media_espacial_media_estacion"]), 1)
            linear_fit_1 = np.poly1d(coef1)
            if (abs(mk.original_test(lista_cldamt[i+2*j]["Media_espacial_media_estacion"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_cldamt[i+2*j]["fecha"], linear_fit_1(
                    np.arange(0, len(lista_cldamt[i+2*j]["fecha"]), 1)), color="k", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_cldamt[i+2*j]["fecha"], linear_fit_1(
                    np.arange(0, len(lista_cldamt[i+2*j]["fecha"]), 1)), color="k", ls=":", lw=0.9)

            coef2 = np.polyfit(np.arange(0, len(lista_bajas[i+2*j]["fecha"]), 1), np.array(
                lista_bajas[i+2*j]["Media_espacial_media_estacion"]), 1)
            linear_fit_2 = np.poly1d(coef2)
            if (abs(mk.original_test(lista_bajas[i+2*j]["Media_espacial_media_estacion"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_bajas[i+2*j]["fecha"], linear_fit_2(np.arange(
                    0, len(lista_bajas[i+2*j]["fecha"]), 1)), color="teal", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_bajas[i+2*j]["fecha"], linear_fit_2(np.arange(
                    0, len(lista_bajas[i+2*j]["fecha"]), 1)), color="teal", ls=":", lw=0.9)

            coef3 = np.polyfit(np.arange(0, len(lista_medias[i+2*j]["fecha"]), 1), np.array(
                lista_medias[i+2*j]["Media_espacial_media_estacion"]), 1)
            linear_fit_3 = np.poly1d(coef3)
            if (abs(mk.original_test(lista_medias[i+2*j]["Media_espacial_media_estacion"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_medias[i+2*j]["fecha"], linear_fit_3(np.arange(
                    0, len(lista_medias[i+2*j]["fecha"]), 1)), color="purple", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_medias[i+2*j]["fecha"], linear_fit_3(np.arange(
                    0, len(lista_medias[i+2*j]["fecha"]), 1)), color="purple", ls=":", lw=0.9)

            coef4 = np.polyfit(np.arange(0, len(lista_altas[i+2*j]["fecha"]), 1), np.array(
                lista_altas[i+2*j]["Media_espacial_media_estacion"]), 1)
            linear_fit_4 = np.poly1d(coef4)
            if (abs(mk.original_test(lista_altas[i+2*j]["Media_espacial_media_estacion"], alpha=0.05)[3]) > 1.96):
                ax[j, i].plot(lista_altas[i+2*j]["fecha"], linear_fit_4(np.arange(
                    0, len(lista_altas[i+2*j]["fecha"]), 1)), color="crimson", ls="-", lw=0.9)
            else:
                ax[j, i].plot(lista_altas[i+2*j]["fecha"], linear_fit_4(np.arange(
                    0, len(lista_altas[i+2*j]["fecha"]), 1)), color="crimson", ls=":", lw=0.9)

            ax[j, i].tick_params(
                axis='x', direction='out', bottom=True, labelrotation=45, labelsize=10, pad=1.5)
            ax[j, i].set_ylim(ymin, ymax)
            ax[j, i].set_xlabel("Fecha", size=10)
            ax[j, i].set_ylabel(variable+" %", size=10)
            ax[j, i].grid()
            ax[j, i].set_title(estacion[j][i])

            ax[j, i].text(1984, 81, "Media (%) \nDesvío (%) \nTendencia (%/dec)",
                          color="k", ha="left", backgroundcolor="white")
            ax[j, i].text(1998, 81, str(np.round(np.nanmean(lista_cldamt[i+2*j]["Media_espacial_media_estacion"]), 1))+" \n"+str(np.round(np.nanstd(
                lista_cldamt[i+2*j]["Media_espacial_media_estacion"]), 2))+" \n"+np.str(np.round(coef1[0]*10, 2)), c="k", backgroundcolor="white")
            ax[j, i].text(2002, 81, str(np.round(np.nanmean(lista_bajas[i+2*j]["Media_espacial_media_estacion"]), 1))+" \n"+str(np.round(np.nanstd(
                lista_bajas[i+2*j]["Media_espacial_media_estacion"]), 2))+" \n"+np.str(np.round(coef2[0]*10, 2)), c="teal", backgroundcolor="white")
            ax[j, i].text(2006, 81, str(np.round(np.nanmean(lista_medias[i+2*j]["Media_espacial_media_estacion"]), 1))+" \n"+str(np.round(np.nanstd(
                lista_medias[i+2*j]["Media_espacial_media_estacion"]), 2))+" \n"+np.str(np.round(coef3[0]*10, 2)), c="purple", backgroundcolor="white")
            ax[j, i].text(2010, 81, str(np.round(np.nanmean(lista_altas[i+2*j]["Media_espacial_media_estacion"]), 1))+" \n"+str(np.round(np.nanstd(
                lista_altas[i+2*j]["Media_espacial_media_estacion"]), 2))+" \n"+np.str(np.round(coef4[0]*10, 2)), c="crimson", backgroundcolor="white")

    fig1.legend(["cldamt", "cldamt bajas", "cldamt medias", "cldamt altas"],
                loc='lower left', ncol=4, bbox_to_anchor=(0.25, 0.00))
    # fig1.tight_layout
    plt.subplots_adjust(left=0.1,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.2,
                        hspace=0.3)
    nombre = variable+"_"+"media_espacial_trimestral_" + \
        region+"_"+"(estaciones)"
    plt.savefig(ruta_salida+nombre, dpi=140)
    plt.show


# %%
ruta_salida = "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_bajas_medias_altas_series/"

serie_trimestral(data_list_cldamt_media_espacial_estacional_sudamerica, data_list_cldamt_bajas_media_espacial_estacional_sudamerica,
                 data_list_cldamt_medias_media_espacial_estacional_sudamerica, data_list_cldamt_altas_media_espacial_estacional_sudamerica, "cldamt", "Sudamérica", 0, 100, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region1, data_list_cldamt_bajas_media_espacial_estacional_region1,
                 data_list_cldamt_medias_media_espacial_estacional_region1, data_list_cldamt_altas_media_espacial_estacional_region1, "cldamt", "Región 1", 0, 100, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_region2, data_list_cldamt_bajas_media_espacial_estacional_region2,
                 data_list_cldamt_medias_media_espacial_estacional_region2, data_list_cldamt_altas_media_espacial_estacional_region2, "cldamt", "Región 2", 0, 100, ruta_salida)
serie_trimestral(data_list_cldamt_media_espacial_estacional_corrientes, data_list_cldamt_bajas_media_espacial_estacional_corrientes,
                 data_list_cldamt_medias_media_espacial_estacional_corrientes, data_list_cldamt_altas_media_espacial_estacional_corrientes, "cldamt", "Corrientes", 0, 100, ruta_salida)

# %%
"""
LISTO CON LAS SERIES, extender a tipos de nube todo el analisis. 
"""
# %%

"""
defino funcion que grafica los campos de una determinada variable en una determinada region
"""


def grafico_campos_nubosidad(paises, provincias, data_list, indice_list, variable, lat_min, lat_max, lon_min, lon_max, unidades_nombre, valor_minimo, valor_maximo, delta_valor, xticks_min, xticks_max, yticks_min, yticks_max, grid, region, ruta_salida, paleta_color):
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
    # carga librerias necesarias
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import cmocean

    # limpio graficos
    plt.close()

    # selecciona variable
    data = data_list[indice_list]
    # selecciona variable y toma el unico valor para cada punto de grilla
    variable_data = data[variable].mean("time", keep_attrs=True)
    # cambio el nombre de la unidad
    variable_data.attrs["units"] = unidades_nombre

    # selecciono region
    lats = variable_data["lat"][:]
    lons = variable_data["lon"][:]
    lat_lims = [lat_min, lat_max]
    lon_lims = [360+lon_min, 360+lon_max]  # lean 360-64 (64 O) 360-31 (31 O)
    lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
    lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
    variable_data_subset = variable_data[lat_inds, lon_inds]

    # extraigo mes
    mes = str(data_list[indice_list]["time"].values[0])[5:7]
    # extraigo anio
    anio = str(data_list[indice_list]["time"].values[0])[0:4]

    # ploteo
    fig1 = plt.figure(figsize=[12, 5], dpi=200)
    ax = fig1.add_subplot(
        111, projection=ccrs.PlateCarree(central_longitude=0))

    if (paleta_color == "rain"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.rain)

    if (paleta_color == "curl"):
        variable_data_subset.plot.contourf(ax=ax,
                                           levels=np.arange(
                                               valor_minimo, valor_maximo, delta_valor),
                                           extend='neither',
                                           transform=ccrs.PlateCarree(),
                                           cbar_kwargs={
                                               'label': variable_data_subset.units},
                                           cmap=cmocean.cm.curl_r)

    ax.add_geometries(provincias, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.5', linewidth=0.7, alpha=0.8)

    ax.add_geometries(paises, crs=ccrs.PlateCarree(), facecolor='none',
                      edgecolor='0.4', alpha=0.8)

    ax.coastlines(color='0.3')

    ax.set_xticklabels(np.arange(xticks_min, xticks_max)[::8])
    plt.xticks(np.arange(xticks_min, xticks_max)[::8])
    ax.set_xlabel("Longitud")

    ax.set_yticklabels(np.arange(yticks_min, yticks_max)[::8])
    plt.yticks(np.arange(yticks_min, yticks_max)[::8])
    ax.set_ylabel("Latitud")

    if (grid == True):
        plt.grid(linestyle="--", alpha=0.3)

    plt.title(variable+" "+region+" "+mes+"/"+anio)
    # plt.tight_layout()
    plt.savefig(ruta_salida+"/"+variable+" "+region+" "+mes+"-"+anio)
    plt.show()


# %%
"""
Grafico los campos de amount de nubosidad (%): cldamt para todos los tiempos y los guardo
"""


# cargo shape con paises

resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'
shpfilename = shapereader.natural_earth(
    resolution, category, name)  # cargo paises de natural_earth
# cargo el shapefile usando geopandas
df = geopandas.read_file(shpfilename)
# leo los paises que voy a necesitar
paises = MultiPolygon([df.loc[df['ADMIN'] == 'Argentina']['geometry'].values[0],
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
                       df.loc[df['ADMIN'] == "Costa Rica"]['geometry'].values[0]])  # los paso a multipolygon para poder graficarlos

# cargo shape con provincias de argentina con datos del IGN
# descargo los datos de aca: https://www.ign.gob.ar/NuestrasActividades/InformacionGeoespacial/CapasSIG "Provincia"
IGN = geopandas.read_file(
    "/home/nadia/Documentos/Doctorado/datos/mapas/provincia/provincia.shp")
provincias = [None]*24
for i in range(0, 24):
    provincias[i] = IGN["geometry"][i]
# paso a multipolygon para poder ponerlo en mapa
provincias = MultiPolygon(provincias)


for i in range(0, cantidad_de_datos):
    grafico_campos_nubosidad(paises, provincias, data_list, i, "cldamt", -39, -16, -64, -31, "%", 0, 101, 5, -60, -31, -
                             35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos", "rain")

# grafico_campos_nubosidad(paises,provincias,data_list,5,"cldamt",-60,15,-90,-30,"%",0,101,5,-85,-30,-55,15,True,"Sudamérica","/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos","rain")

# %%
"""

Calculo anomalias climaticas mensuales de variable. Tomo los anios enteros: 1984 a 2016 (32 anios) (DUDA seria hasta 2014?)

1) Calculo media climatica de variable mensual para cada mes, es decir por ejemplo para Enero seria el promedio de la cantidad punto a punto de todos los Eneros del periodo climatologico.

2) A cada mes de cada anio punto a punto le resto la media de ese mes calculada en 1

"""


# %%


# %%
# hago un array 3d en donde cada capa es un campo mensual y calculo la media
def media_mensual(data_list, variable, mes):
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
    n = 1
    arr_3D = np.empty((1, 180, 360))
    for i in range(0, len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7] == mes):
            variable_data = [data_list[i][variable].mean(
                "time", keep_attrs=True).values]
            arr_3D = np.concatenate([arr_3D, variable_data])
            n = n+1
            arr_3D = np.reshape(arr_3D, (n, 180, 360))
    arr_3D = arr_3D[1:np.shape(arr_3D)[0], :, :]
    media_mensual = np.mean(arr_3D, axis=0)
    return(media_mensual)

# calculo la anomalia mensual restandole a cada mes la media de ese mes y lo agrego al xarray correspondiente de la lista data_list con nombre "media_climatologica_"+variable


def anomalia_mensual(data_list, variable, mes):
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
    for i in range(0, len(data_list)):
        if (str(data_list[i]["time"].values[0])[5:7] == mes):
            variable_data = [data_list[i][variable].mean(
                "time", keep_attrs=True).values]
            anom = variable_data-media_mensual(data_list, variable, mes)
            anom_dataarray = xr.DataArray(
                data=anom, dims=["time", "lat", "lon"])
            data_list[i] = data_list[i].assign(variable_anom=anom_dataarray)
            data_list[i] = data_list[i].rename(
                {"variable_anom": "anomalia_mensual_"+variable})
    return(data_list)


data_list_modificado = data_list.copy()
meses = ["01", "02", "03", "04", "05", "06",
         "07", "08", "09", "10", "11", "12"]
for i in range(0, 12):
    data_list_modificado = anomalia_mensual(
        data_list_modificado, "cldamt", meses[i])


for i in range(0, cantidad_de_datos):
    grafico_campos_nubosidad(paises, provincias, data_list_modificado, i, "anomalia_mensual_cldamt", -39, -16, -64, -31, "%", -60, 65, 5, -60, -31, -
                             35, -18, True, "Región 1", "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_anomalias_mensuales_campos", "curl")

# %%

"""
calculo la media mensual de una determinada variable en determinada region 
grafico una serie temporal de todo el periodo, y por mes
"""
# defino funcion que calcula media mensual de una determinada variable en determinada region


def media_espacial(data_list, indice_list, variable, lat_min, lat_max, lon_min, lon_max):
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

    data = data_list[indice_list]
    # selecciona variable y toma el unico valor para cada punto de grilla
    variable_data = data[variable].mean("time", keep_attrs=True)

    # selecciono region
    lats = variable_data["lat"][:]
    lons = variable_data["lon"][:]
    lat_lims = [lat_min, lat_max]
    lon_lims = [360+lon_min, 360+lon_max]  # lean 360-64 (64 O) 360-31 (31 O)
    lat_inds = np.where((lats > lat_lims[0]) & (lats < lat_lims[1]))[0]
    lon_inds = np.where((lons > lon_lims[0]) & (lons < lon_lims[1]))[0]
    variable_data_subset = variable_data[lat_inds, lon_inds]

    # calculo media espacial del subset de la variable data
    media_espacial = np.mean(variable_data_subset.values, axis=(0, 1))
    return(media_espacial)

# armo funcion que devuelve data frame donde primera columna sea la fecha y la segunda columna sea la media espacial


def media_espacial_df(data_list, variable, lat_min, lat_max, lon_min, lon_max):
    media_espacial_df = pd.DataFrame(
        columns=["fecha", "Media_espacial_"+variable])
    for i in range(0, len(data_list)):
        # extraigo mes
        mes = str(data_list[i]["time"].values[0])[5:7]
        # extraigo anio
        anio = str(data_list[i]["time"].values[0])[0:4]
        # ¢alculo media
        media_espacial_i = media_espacial(
            data_list, i, variable, lat_min, lat_max, lon_min, lon_max)
        media_espacial_df = media_espacial_df.append(
            {"fecha": mes+"-"+anio, "Media_espacial_"+variable: media_espacial_i}, ignore_index=True)

    return(media_espacial_df)


cldamt_media_espacial_df = media_espacial_df(
    data_list, "cldamt", -60, -31, -35, -18)


fig1, ax = plt.subplots(figsize=[12, 6], dpi=200)
plt.plot(cldamt_media_espacial_df["fecha"],
         cldamt_media_espacial_df["Media_espacial_cldamt"])
major_ticksx = np.arange(6, len(cldamt_media_espacial_df["fecha"]), len(
    cldamt_media_espacial_df["fecha"])/34)
ax.set_xticks(major_ticksx)
ax.grid(alpha=0.3)
ax.tick_params(axis='x', direction='out', bottom=True,
               labelrotation=45, labelsize=10, pad=1.5)
ax.set_xlabel("Fecha", size=10)
ax.set_ylabel("cldamt media %", size=10)
plt.title("Media espacial mensual cldmt region 1")
plt.savefig(
    "/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos/serie", dpi=140)
