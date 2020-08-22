import itertools
from collections import Counter
import sys

from lsst.sims.maf.metrics import BaseMetric
import lsst.sims.maf.metrics as metrics
import healpy as hp
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta

# sys.path.append("/home/idies/workspace/Storage/aagonzalez6/persistent/LSST_OpSim/Scripts_NBs/Essentials/")
from .Lookup_table_read import Ntot_lookup


# si se quiere cambiar el rango z de integracion para Ntot hay que modificar el archivo
# lookup que esta dentro de Lookup_table_read.py
def LogHist(data,nlogbins=50,exponentes_limites = "auto"):
    """ data = 1d lista
        nlogbins = int, cantidad de bines en los que agrupar los datos
        exponentes_limites = [lim_inf,lim_sup] o "auto" 
        
        Entrega array([valores y]), array([bins])"""
    

    if exponentes_limites == "auto":# ______________________________"auto"_________________________
    
        if max(data)==0:
            cota_exponente_sup = 1
        else:
            cota_exponente_sup = int(np.log10(max(data)))+1

        if min(data)<=0:
            binsitos = np.logspace(0,cota_exponente_sup,nlogbins) 

        else:
            cota_exponente_inf = int(np.log10(min(data)))
            binsitos = np.logspace(cota_exponente_inf,cota_exponente_sup,nlogbins) 
    else: #_______________________________________________limites exponentes dados _________________
        binsitos = np.logspace(exponentes_limites[0],exponentes_limites[1],nlogbins)
        
        
    return np.histogram(data,binsitos) # eso es ([valBin1,valBin2,...,valBinN],[])
    
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################
class NpairsMetric(BaseMetric): # el backup antes de 5 de agosto esta en sciserver
    def __init__(self, nightCol='night', metricName='NpairsMetric', units = 'Cantidad'):#, metricDtype="object"):
#    def __init__(self, mjdTimeCol='observationStartMJD', metricName='NpairsMetric', units = 'Cantidad'):#, metricDtype="object"): #Version segundos
        self.nightCol = nightCol 
        super(NpairsMetric, self).__init__(col=nightCol, metricName=metricName, units = units, metricDtype="object")
#        super(NpairsMetric, self).__init__(col=mjdTimeCol, metricName=metricName, units = units, metricDtype="object") #Version segundos
#	 self.mjdTimeCol = mjdTimeCol #Version segundos

    def run(self, dataSlice, slicePoint):
        array_noches = dataSlice[self.nightCol]
#        array_mjdTime = dataSlice[self.mjdTimeCol] #Version segundos
        """Metodo 1"""
#         pares = itertools.combinations(array_noches,2)

#         dictDiferencias = Counter([np.abs(x[1]-x[0]) for x in pares]) # cantidad de diferencia temporal en dias
#         dict_new = {k:v for k,v in dictDiferencias.items() if  k == 7}
#         return dictDiferencias
        """Metodo 2""" # entrega lista con diferecias nomas, dsp hay que contar cuantas dif iguales (que estan en dias) hay/ 
        
        lista_distancias_pares = pdist(array_noches[:, None], 'cityblock').astype("int")
#        lista_distancias_pares = pdist(array_mjdTime[:, None], 'cityblock').astype("float") #Version segundos      
#        lista_distancias_pares = np.array([TimeDelta(x, format="jd").sec for x in lista_distancias_pares]) #Version segundos
        
        if len(lista_distancias_pares)>0:
            return LogHist(lista_distancias_pares,exponentes_limites = [0,4])
        else:
            return None
    
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################
class NtotMetricV2(BaseMetric):
    def __init__(self, mod, f1f2diff = 2,fivesigmacol='fiveSigmaDepth', metricName='NtotMetricV2', units = 'Cantidad'):
        ver_int = "V5_" #### si se ocupa version 4 ocupar en blanco:  ""
        self.lookup_tables = {"u":"{}u_mod{}_LookupT_extension.pkl".format(ver_int,mod),
                            "g":"{}g_mod{}_LookupT_extension.pkl".format(ver_int,mod),
                            "r":"{}r_mod{}_LookupT_extension.pkl".format(ver_int,mod),
                            "i":"{}i_mod{}_LookupT_extension.pkl".format(ver_int,mod),
                            "z":"{}z_mod{}_LookupT_extension.pkl".format(ver_int,mod),
                            "y":"{}y_mod{}_LookupT_extension.pkl".format(ver_int,mod)}#debe ser de la forma de {"u_lookuptable":}
        
#         self.f1f2 = f1f2 # tupla (filter1,filter2) e.g.("z","y") # DEBE estar en ORDEN de mas azul a mas rojo
        self.coaddedm5 = metrics.Coaddm5Metric(m5Col = fivesigmacol) # este m5col no era necesario (viene por defecto)
        self.Nsigma = 10 #Nsigma. En este caso 10 sigma
        self.f1f2diff = f1f2diff

        super(NtotMetricV2, self).__init__(col=[fivesigmacol,"filter"], metricName=metricName, units = units)
    
    def dict_f1f2(self, f1, f2, dataSlice_f1band, dataSlice_f2band, area):
        lookup_table_f2 = self.lookup_tables[f2] # el nombre del archivo que debe ubicarse en .../Essentials/LookupTables/ 
        
        if len(dataSlice_f2band) == 0:# si pasa de este if entonces dataSlice_f2band contiene elementos
            return {f2+"5_Ntot":0,
                    f2+"10_Ntot":0,
                    f1+f2+"5_Ntot":0,
                    f1+f2+"10_Ntot":0}
        f2_coadd5sigma_this_Hpix = self.coaddedm5.run(dataSlice_f2band) 
        f2_coaddNsigma_this_Hpix = f2_coadd5sigma_this_Hpix - np.log10(self.Nsigma/5)
        
        if len(dataSlice_f1band) == 0 : # si pasa de este if entonces f2band y f1band contienen elementos
            return {f2+"5_Ntot":Ntot_lookup(f2_coadd5sigma_this_Hpix,area,lookup_table_f2),
                    f2+"10_Ntot":Ntot_lookup(f2_coaddNsigma_this_Hpix,area,lookup_table_f2),
                    f1+f2+"5_Ntot":0,
                    f1+f2+"10_Ntot":0} # no hay elementos f1band para aplicar f1f2_Ntot
        f1_coadd5sigma_this_Hpix = self.coaddedm5.run(dataSlice_f2band)
        f1_coaddNsigma_this_Hpix = f1_coadd5sigma_this_Hpix - np.log10(self.Nsigma/5)

        if f1_coaddNsigma_this_Hpix - f2_coaddNsigma_this_Hpix >= self.f1f2diff:
            return {f2+"5_Ntot":Ntot_lookup(f2_coadd5sigma_this_Hpix,area,lookup_table_f2),
                    f2+"10_Ntot":Ntot_lookup(f2_coaddNsigma_this_Hpix,area,lookup_table_f2),
                    f1+f2+"5_Ntot":Ntot_lookup(f2_coadd5sigma_this_Hpix,area,lookup_table_f2),
                    f1+f2+"10_Ntot":Ntot_lookup(f2_coaddNsigma_this_Hpix,area,lookup_table_f2)} ## si se cambia pa dar opcion al lookup hay que modificar aqui tmb --ya lo hice
        else: # en este caso hay un rango de magnitudes y detectables que al estar menos alejadas de z no podrian ser confirmadas, por tanto solo se toman en cuentas las magnitudes confirmables ie y = z-2 para abajo
            return {f2+"5_Ntot":Ntot_lookup(f2_coadd5sigma_this_Hpix,area,lookup_table_f2),
                    f2+"10_Ntot":Ntot_lookup(f2_coaddNsigma_this_Hpix,area,lookup_table_f2),
                    f1+f2+"5_Ntot":Ntot_lookup(f1_coadd5sigma_this_Hpix - 2 ,area,lookup_table_f2),
                    f1+f2+"10_Ntot":Ntot_lookup(f1_coaddNsigma_this_Hpix - 2 ,area,lookup_table_f2)}

    def run(self, dataSlice, slicePoint): ######## DATASLICE ES UNA STRUCUTURED ARRAY !! ES "COMO" UN DICT PERO POR COLUMNA
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)
# dataSlice e.g.:        
#         fiveSigmaDepth  |  filter  |...| otros
#         ----------------------------...-------
#                24       |    z     |
#                25       |    z     |
#                23.3     |    y     |             
#                22.1     |    z     |                
#                         .          .
#                         .          .
#                         .          .
        
        mask_u = (dataSlice["filter"] == "u")
        mask_g = (dataSlice["filter"] == "g")
        mask_r = (dataSlice["filter"] == "r") 
        mask_i = (dataSlice["filter"] == "i") 
        mask_z = (dataSlice["filter"] == "z")    
        mask_y = (dataSlice["filter"] == "y")

        dataSlice_u = dataSlice[mask_u]
        dataSlice_g = dataSlice[mask_g]
        dataSlice_r = dataSlice[mask_r]
        dataSlice_i = dataSlice[mask_i]
        dataSlice_z = dataSlice[mask_z]
        dataSlice_y = dataSlice[mask_y]
        
        lista_diccionarios = []
        dict_u = self.dict_f1f2("x", "u", [],          dataSlice_u, area)
        lista_diccionarios.append(dict_u)
        dict_g = self.dict_f1f2("u", "g", dataSlice_u, dataSlice_g, area)
        lista_diccionarios.append(dict_g)
        dict_r = self.dict_f1f2("g", "r", dataSlice_g, dataSlice_r, area)
        lista_diccionarios.append(dict_r)
        dict_i = self.dict_f1f2("r", "i", dataSlice_r, dataSlice_i, area)
        lista_diccionarios.append(dict_i)
        dict_z = self.dict_f1f2("i", "z", dataSlice_i, dataSlice_z, area)
        lista_diccionarios.append(dict_z)
        dict_y = self.dict_f1f2("z", "y", dataSlice_z, dataSlice_y, area)
        lista_diccionarios.append(dict_y)
        dict_area = {"Area":area }
        lista_diccionarios.append(dict_area)

        # dict maestro va a tener 4*6 (ntots) +1 (area) elementos
        # lista_diccionarios = [dict_u,.. dict_y, dict_area]

        return lista_diccionarios # este es lo que reducefunction toman como metricValue
        

    def reduceu_o5_Ntot(self, metricValue):return metricValue[0]["u5_Ntot"]
    def reduceu_o10_Ntot(self, metricValue):return metricValue[0]["u10_Ntot"]
    def reducexu_o5_Ntot(self, metricValue):return metricValue[0]["xu5_Ntot"]
    def reducexu_o10_Ntot(self, metricValue):return metricValue[0]["xu10_Ntot"]
    # ====================================================
        
    def reduceg_o5_Ntot(self, metricValue):return metricValue[1]["g5_Ntot"]
    def reduceg_o10_Ntot(self, metricValue):return metricValue[1]["g10_Ntot"]
    def reduceug_o5_Ntot(self, metricValue):return metricValue[1]["ug5_Ntot"]
    def reduceug_o10_Ntot(self, metricValue):return metricValue[1]["ug10_Ntot"]
    # ====================================================
    
    def reducer_o5_Ntot(self, metricValue):return metricValue[2]["r5_Ntot"]
    def reducer_o10_Ntot(self, metricValue):return metricValue[2]["r10_Ntot"]
    def reducegr_o5_Ntot(self, metricValue):return metricValue[2]["gr5_Ntot"]
    def reducegr_o10_Ntot(self, metricValue):return metricValue[2]["gr10_Ntot"]
    # ====================================================
    
    def reducei_o5_Ntot(self, metricValue):return metricValue[3]["i5_Ntot"]
    def reducei_o10_Ntot(self, metricValue):return metricValue[3]["i10_Ntot"]
    def reduceri_o5_Ntot(self, metricValue):return metricValue[3]["ri5_Ntot"]
    def reduceri_o10_Ntot(self, metricValue):return metricValue[3]["ri10_Ntot"]
    # ====================================================
    
    def reducez_o5_Ntot(self, metricValue):return metricValue[4]["z5_Ntot"]
    def reducez_o10_Ntot(self, metricValue):return metricValue[4]["z10_Ntot"]
    def reduceiz_o5_Ntot(self, metricValue):return metricValue[4]["iz5_Ntot"]
    def reduceiz_o10_Ntot(self, metricValue):return metricValue[4]["iz10_Ntot"]
    # ====================================================
    
    def reducey_o5_Ntot(self, metricValue):return metricValue[5]["y5_Ntot"]
    def reducey_o10_Ntot(self, metricValue):return metricValue[5]["y10_Ntot"]
    def reducezy_o5_Ntot(self, metricValue):return metricValue[5]["zy5_Ntot"]
    def reducezy_o10_Ntot(self, metricValue):return metricValue[5]["zy10_Ntot"]
    # ====================================================
    def reducedict_area(self, metricValue):
         return metricValue[6]["Area"]         

        
        
        
