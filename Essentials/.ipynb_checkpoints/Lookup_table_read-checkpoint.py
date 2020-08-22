import pickle
import numpy as np
from bisect import bisect_left
import os

actual_dir = os.path.dirname(os.path.abspath(__file__))

# ========================== Nearest key function ===================================================
def nearest_key(valor,numDict): 
	"""
	valor = float
	numDic = diccionario numerico (con keys = float)

	Entrega la key más cercana a valor. Si hay dos keys igual de cerca, se elige la key más pequeña.
	"""

	lista_keys = list(numDict.keys())
	i = bisect_left(lista_keys,valor) # todos los elementos de la lista con indice <i son <valor
	posibles_3 = [lista_keys[i-1], lista_keys[i], lista_keys[i+1]]
	return min(posibles_3, key = lambda x: abs(x - valor)) # si i-1 e i+1 tienen igual diferencia, se elige el i-1
	
# =============================================================================
def Ntot_lookup(my,omega,archivo):
    """
    my = limiting magnitude Y band
    omega = area de obs (radians)
    archivo = str, nombre del archivo e.g. hola.pkl
    
    """
    # ========================= Importar diccionario =============================
    a_file = open(actual_dir+"/LookupTables/{}".format(archivo), "rb") #z=5.7 a 6.5
    lookup_dict = pickle.load(a_file)
    a_file.close()
    # ============================================================================
    
    return omega*lookup_dict[nearest_key(my,lookup_dict)]


