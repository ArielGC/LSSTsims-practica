import itertools
from collections import Counter
import sys

from lsst.sims.maf.metrics import BaseMetric
import lsst.sims.maf.metrics as metrics
#from Integral import Ntot
import healpy as hp
import numpy as np
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt

sys.path.append("/home/idies/workspace/Storage/aagonzalez6/persistent/LSST_OpSim/Scripts_NBs/Essentials/")
from Lookup_table_read import Ntot_lookup

#from quasarlf.pubtools.utilities import return_qlf_in_band

#print(return_qlf_in_band(2, -5, model='A'))
#def hola():
#    return return_qlf_in_band(2, -5, model='A')


# >>> z = ['blue', 'red', 'blue', 'yellow', 'blue', 'red']
# >>> Counter(z)
# Counter({'blue': 3, 'red': 2, 'yellow': 1})

# si se quiere cambiar el rango z de integracion para Ntot hay que modificar el archivo
# lookup que esta dentro de Lookup_table_read.py
def LogHist(data,nlogbins=50,exponentes_limites = "auto"):
    """ data = 1d lista
        nlogbins = int, cantidad de bines en los que agrupar los datos
        exponentes_limites = [lim_inf,lim_sup] o "auto" 
        
        Entrega array([valores y]), array([bins])"""
    
    #potencia de 10 mas cercana
#     try:cota_exponente_sup = int(np.log10(max(data)))+1
#     except: raise ValueError(data)
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

class NQsoMetric5m(BaseMetric):
    """Calculate the coadded m5 value at this gridpoint."""
    def __init__(self,lookup_table, m5col='fiveSigmaDepth', metricName='NQsoMetric5m', units = 'Cantidad'):
        """Instantiate metric.
        m5col = the column name of the individual visit m5 data."""
        self.lookup_table = lookup_table
        self.m5col = m5col
        super(NQsoMetric5m, self).__init__(col=m5col, metricName=metricName, units = units)

    def run(self, dataSlice, slicePoint):
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)
        coadd5m_this_Hpix = 1.25 * np.log10(np.sum(10.**(.8*dataSlice[self.m5col])))
       
        return Ntot_lookup(coadd5m_this_Hpix,area,self.lookup_table)
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################

class NQsoMetric10m(BaseMetric):
    """Calculate the coadded m5 value at this gridpoint."""
    def __init__(self,lookup_table, m5col='fiveSigmaDepth', metricName='NQsoMetric10m', units = 'Cantidad'):
        """Instantiate metric.
        m5col = the column name of the individual visit m5 data."""
        self.lookup_table=lookup_table
        self.m5col = m5col
        super(NQsoMetric10m, self).__init__(col=m5col, metricName=metricName, units = units)

    def run(self, dataSlice, slicePoint):
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)
        coadd5m_this_Hpix = 1.25 * np.log10(np.sum(10.**(.8*dataSlice[self.m5col])))
        coadd10m_this_Hpix = coadd5m_this_Hpix - np.log10(2)
       
        return Ntot_lookup(coadd10m_this_Hpix,area,self.lookup_table)

class AreaTot(BaseMetric):
    """Calculate the coadded m5 value at this gridpoint."""
    def __init__(self, m5col='fiveSigmaDepth', metricName='Area', units = 'Cantidad'):
        """Instantiate metric.
        m5col = the column name of the individual visit m5 data."""
        self.m5col = m5col
        super(AreaTot, self).__init__(col=m5col, metricName=metricName, units = units)

    def run(self, dataSlice, slicePoint):
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)       
        return area
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################

class NQsoMetric10mv2(BaseMetric):
    def __init__(self,lookup_table, colsita='fiveSigmaDepth', metricName='NQsoMetric10mv2', units = 'Cantidad'):
        """Instantiate metric.
        m5col = the column name of the individual visit m5 data."""
        self.lookup_table = lookup_table
        self.coaddedm5 = metrics.Coaddm5Metric(m5Col = colsita) # este m5col no era necesario (viene por defecto)
        # you must pass the columns you 
        # need, and the kwargs to the call of super.
        super(NQsoMetric10mv2, self).__init__(col=colsita, metricName=metricName, units = units)

    def run(self, dataSlice, slicePoint):
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)
        coadd5m_this_Hpix = self.coaddedm5.run(dataSlice)
        coadd10m_this_Hpix = coadd5m_this_Hpix - np.log10(2)
       
        return {"Ntot10m":Ntot_lookup(coadd10m_this_Hpix,area,lookup_table),"Area":area,"DataSlice10mv2":len(dataSlice)}
    
    def reduceNtot321(self, metricValue):
        return metricValue["Ntot10m"]    
    
    def reduceArea(self, metricValue):
        return metricValue["Area"]
    def reduceDataSlice10mv2(self, metricValue):
        return metricValue["DataSlice10mv2"]
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################
class NpairsMetric(BaseMetric):
    def __init__(self, nightCol='night', metricName='NpairsMetric', units = 'Cantidad'):#, metricDtype="object"):
        self.nightCol = nightCol
        super(NpairsMetric, self).__init__(col=nightCol, metricName=metricName, units = units, metricDtype="object")

    def run(self, dataSlice, slicePoint):
        array_noches = dataSlice[self.nightCol]
        """Metodo 1"""
#         pares = itertools.combinations(array_noches,2)

#         dictDiferencias = Counter([np.abs(x[1]-x[0]) for x in pares]) # cantidad de diferencia temporal en dias
#         dict_new = {k:v for k,v in dictDiferencias.items() if  k == 7}
#         return dictDiferencias
        """Metodo 2""" # entrega lista con diferecias nomas, dsp hay que contar cuantas dif iguales (que estan en dias) hay/ 
              
        lista_distancias_pares = pdist(array_noches[:, None], 'cityblock').astype("int")      
        #dictDiferencias = Counter(lista_distancias_pares)
        #dict_new = {k:v for k,v in dictDiferencias.items() if }
        if len(lista_distancias_pares)>0:
            return LogHist(lista_distancias_pares,exponentes_limites = [0,4])#len(lista_distancias_pares) #dict_new
        else:
            return None
    
    # si da error con esto,probar con metricdtype cambiar a integer, si da nuevamente error
    # es por el dict Diferencias. Si no, era por las combinaciones
    
# ####################################################################################################
# ####################################################################################################
# ****************************************************************************************************
# ####################################################################################################
# ####################################################################################################

class NtotMetricZY(BaseMetric):
    def __init__(self,lookup_table, zydiff = 2,fivesigmacol='fiveSigmaDepth', metricName='NtotMetricZY', units = 'Cantidad'):
        self.lookup_table = lookup_table
        self.coaddedm5 = metrics.Coaddm5Metric(m5Col = fivesigmacol) # este m5col no era necesario (viene por defecto)
        self.Nsigma = 10 #Nsigma
        self.zydiff = zydiff
        super(NtotMetricZY, self).__init__(col=[fivesigmacol,"filter"], metricName=metricName, units = units)
    

    def run(self, dataSlice, slicePoint): ######## DATASLICE ES UNA STRUCUTURED ARRAY !! ES "COMO" UN DICT PERO POR COLUMNA
        nside = slicePoint["nside"]
        area = hp.nside2pixarea(nside)
# dataSlice:        
#         fiveSigmaDepth  |  filter  |...| otros
#         ----------------------------...-------
#                24       |    z     |
#                25       |    z     |
#                23.3     |    y     |             
#                22.1     |    z     |                
#                         .          .
#                         .          .
#                         .          .
                        
        maskz = (dataSlice["filter"] == "z")
        masky = (dataSlice["filter"] == "y")
        dataSlice_zband = dataSlice[maskz]
        dataSlice_yband = dataSlice[masky] #
        
        if len(dataSlice_yband) == 0:# si pasa de este if entonces dataSlice_yband contiene elementos
            return {"y5_Ntot":0,"y10_Ntot":0,"zy5_Ntot":0,"zy10_Ntot":0,"Area":area}
        y_coadd5sigma_this_Hpix = self.coaddedm5.run(dataSlice_yband) 
        y_coaddNsigma_this_Hpix = y_coadd5sigma_this_Hpix - np.log10(self.Nsigma/5)
        
        if len(dataSlice_zband) == 0: # si pasa de este if entonces yband y zband contienen elementos
            return {"y5_Ntot":Ntot_lookup(y_coadd5sigma_this_Hpix,area,self.lookup_table),
                    "y10_Ntot":Ntot_lookup(y_coaddNsigma_this_Hpix,area,self.lookup_table),
                    "zy5_Ntot":0,
                    "zy10_Ntot":0,
                    "Area":area} # no hay elementos zband para aplicar zy_Ntot
        z_coadd5sigma_this_Hpix = self.coaddedm5.run(dataSlice_yband)
        z_coaddNsigma_this_Hpix = z_coadd5sigma_this_Hpix - np.log10(self.Nsigma/5)

        if z_coaddNsigma_this_Hpix - y_coaddNsigma_this_Hpix >= self.zydiff:
            return {"y5_Ntot":Ntot_lookup(y_coadd5sigma_this_Hpix,area,self.lookup_table),
                    "y10_Ntot":Ntot_lookup(y_coaddNsigma_this_Hpix,area,self.lookup_table),
                    "zy5_Ntot":Ntot_lookup(y_coadd5sigma_this_Hpix,area,self.lookup_table),
                    "zy10_Ntot":Ntot_lookup(y_coaddNsigma_this_Hpix,area,self.lookup_table),
                    "Area":area} ## si se cambia pa dar opcion al lookup hay que modificar aqui tmb --ya lo hice
        else: # en este caso hay un rango de magnitudes y detectables que al estar menos alejadas de z no podrian ser confirmadas, por tanto solo se toman en cuentas las magnitudes confirmables ie y = z-2 para abajo
            return {"y5_Ntot":Ntot_lookup(y_coadd5sigma_this_Hpix,area,self.lookup_table),
                    "y10_Ntot":Ntot_lookup(y_coaddNsigma_this_Hpix,area,self.lookup_table),
                    "z5_Ntot":Ntot_lookup(z_coadd5sigma_this_Hpix - 2 ,area,self.lookup_table),
                    "zy10_Ntot":Ntot_lookup(z_coaddNsigma_this_Hpix - 2 ,area,self.lookup_table),
                    "Area":area}            

    def reducey5_Ntot(self, metricValue):
        return metricValue["y5_Ntot"]   
    def reducey10_Ntot(self, metricValue):
        return metricValue["y10_Ntot"]  
    
    def reducezy5_Ntot(self, metricValue):
        return metricValue["zy5_Ntot"]  
    def reducezy10_Ntot(self, metricValue):
        return metricValue["zy10_Ntot"] 
    
    def reduceArea(self, metricValue):
        return metricValue["Area"]
    
#######################################################################################################################
#######################################################################################################################

class NtotMetricV2(BaseMetric):
    def __init__(self, mod, f1f2diff = 2,fivesigmacol='fiveSigmaDepth', metricName='NtotMetricV2', units = 'Cantidad'):
        self.lookup_tables = {"u":"u_mod{}_LookupT_extension.pkl".format(mod),
                            "g":"g_mod{}_LookupT_extension.pkl".format(mod),
                            "r":"r_mod{}_LookupT_extension.pkl".format(mod),
                            "i":"i_mod{}_LookupT_extension.pkl".format(mod),
                            "z":"z_mod{}_LookupT_extension.pkl".format(mod),
                            "y":"y_mod{}_LookupT_extension.pkl".format(mod)}#debe ser de la forma de {"u_lookuptable":}
        
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

        
        
        
        
#     def reducedict_xu(self, metricValue):
#         return metricValue[0]         

#     def reducedict_ug(self, metricValue):
#         return metricValue[1]     

#     def reducedict_gr(self, metricValue):
#         return metricValue[2]     

#     def reducedict_ri(self, metricValue):
#         return metricValue[3]     

#     def reducedict_iz(self, metricValue):
#         return metricValue[4]  

#     def reducedict_zy(self, metricValue):
#         return metricValue[5]       

#     def reducedict_area(self, metricValue):
#         return metricValue[6]["Area"]         



