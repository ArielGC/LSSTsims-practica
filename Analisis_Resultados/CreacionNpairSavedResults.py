import os
import inspect
import glob
import sys
import pickle
# sys.path.append("/home/idies/workspace/Storage/aagonzalez6/persistent/LSST_OpSim/Scripts_NBs/Essentials")

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
# import convenience functions
from Essentials.opsimUtils import *
from Essentials.funcionesAnalisis import HiddenPrints
from Essentials.FMNtot import NpairsMetric

def SaveFile(datos,nombre_guardado):
    dir_guardado = "/data/agonzalez/NpairFIXsamefilt/NpairSavedResults"
    a_file = open(dir_guardado + "/{}".format(nombre_guardado), "wb")
    pickle.dump(datos, a_file)
    a_file.close()
    print("Done:", nombre_guardado)
    


# In[22]:


def analisis_en_masa(FBS, filtro, modo = "save"):

    # ==========================================================        
    filtros_considerados = [filtro]
    pruebas = "" #escribir PRUEBAS si si, si no, dejar en blanco
    # ==========================================================
    
    resultDbPath = '/data/agonzalez/NpairFIXsamefilt/output_FBS_{0}{2}/Npairs_{1}/'.format(FBS,"".join(filtros_considerados),pruebas)
    metricDataPath = '/data/agonzalez/NpairFIXsamefilt/output_FBS_{0}{2}/Npairs_{1}/MetricData/'.format(FBS,"".join(filtros_considerados),pruebas)
    # get a dictionary of resultDb from given directory
    resultDbs = getResultsDbs(resultDbPath)
    # print(resultDbs)
    # the following line will be useful if you did not run MAF on all 75 opsims
    runNames = list(resultDbs.keys())
    # retrieve metricBundles for each opsim run and store them in a dictionary
    bundleDicts = {} # metric bundles de todas las opsim : {runName1:bundledict1,...}
    
    registro = []
    for runName in runNames:
        with HiddenPrints():
            bd_opsim = bundleDictFromDisk(resultDbs[runName], runName, metricDataPath)
        bundleDicts[runName] = bd_opsim

        reg_opsim = ["---","---",runName] #(DDF,WFD,runName)
        for x in list(bd_opsim.items()):
#             print(x[1].constraint)
#             print("filter = '{}' and proposalId > 1".format(filtro))
            if x[1].constraint == "filter = '{}' and note LIKE '%DD%'".format(filtro):  # edit luego de finalizar todo--- acordarse de lo que pasaba con constraint "note LIKE '%DD%'"
                bundle_DDF = x[1]
                data_DDF = bundle_DDF.metricValues.data # esto es una lista de np hists asociados a
                reg_opsim[0] = "DDF"
                if modo == "save":
                    SaveFile(data_DDF,"Npair_data_FBS{0}_{1}_DDF@{2}.pkl".format(FBS, filtro, runName))
            elif x[1].constraint == "filter = '{}' and note NOT LIKE '%DD%'".format(filtro): # la constraint es esta , pero solo
                #porque no fue necesario ponerla completa ya que se trabajo sobre el healpixslicer para aislar
                #la WFD area
                bundle_WFD = x[1] # en vola fuente error porquedsp pongo una lista que usa estos valores sin estar 100% seguro si existiran
                data_WFD = bundle_WFD.metricValues.data
                reg_opsim[1] = "WFD"
                if modo == "save":
                    SaveFile(data_WFD,"Npair_data_FBS{0}_{1}_WFD@{2}.pkl".format(FBS, filtro, runName))
        registro.append(reg_opsim)
        print("---------- Done FBS{}, filtro {}, opsim {}".format(FBS,filtro,runName))
    return registro


# In[23]:


registro_general = []
for FBS in ["1.4.1","1.5"]:
    registro_FBS = []
    for filtro in ["u","g","r","i","z","y"]:
        registro_FBS += analisis_en_masa(FBS, filtro, modo = "save" )
        print("========== done filtro {} de {}".format(filtro,FBS))
    print("========================== done  FBS {}".format(FBS))
    registro_FBS += ["##################################################################################"]
    registro_general += registro_FBS
print("==========================================DONE ALL================================")
print ("__________________________________________________________________________________________")
print ("|||||||||||||||||||||||||||||| PRINTEANDO REGISTRO GENERAL: ||||||||||||||||||||||||||||||")
for x in registro_general: print(x)

