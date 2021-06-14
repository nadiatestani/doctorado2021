"""
Codigo que carga imagenes y arma gif
"""

def gif(ruta,nombre_comun_carga,nombre_comun_salida,nombre_comun_diferencia):
    import glob
    from PIL import Image
    
    fp_in = ruta+nombre_comun_carga+nombre_comun_diferencia+"*.png"
    fp_out=ruta+nombre_comun_salida+nombre_comun_diferencia+".gif"
    
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(fp=fp_out, format='GIF', append_images=imgs,save_all=True, duration=700, loop=0)


#%% Armo GIF de cldamt y de anomalias

ruta="/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_campos/"
nombre_comun_carga="cldamt Región 1 "
nombre_comun_salida="cldamt Región 1 "
nombre_comun_diferencia_vector=["01", "02", "03", "04", "05","06","07","08","09","10","11","12"]

for i in range(0,12):
    gif(ruta,nombre_comun_carga,nombre_comun_salida,nombre_comun_diferencia_vector[i])



ruta="/home/nadia/Documentos/Doctorado/resultados/resultados2021/nubosidad/cldamt_anomalias_mensuales_campos/"
nombre_comun_carga="anomalia_mensual_cldamt Región 1 "
nombre_comun_salida="cldamt_mensual_cldamt Región 1 "
nombre_comun_diferencia_vector=["01", "02", "03", "04", "05","06","07","08","09","10","11","12"]

for i in range(0,12):
    gif(ruta,nombre_comun_carga,nombre_comun_salida,nombre_comun_diferencia_vector[i])
