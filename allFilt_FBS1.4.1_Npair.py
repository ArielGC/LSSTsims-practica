#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os
import inspect
import glob
import sys
# sys.path.append("/home/idies/workspace/Storage/aagonzalez6/persistent/LSST_OpSim/Scripts_NBs/Essentials") # sciserver line

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
#%matplotlib inline #con nbconvert esto daba problema

from Essentials.FMNtot import NpairsMetric

# In[2]:


# For the WFD footprint definitions
from lsst.sims.featureScheduler import utils as schedUtils
print(schedUtils.__file__)
# To calculate metrics with MAF
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as mb
import lsst.sims.maf.plots as plots
import lsst.sims.maf.batches as batches
# import convenience functions
from Essentials.opsimUtils import *
# ______________________________________________________________________________________________________

# In[3]:


lista_filtros = ["u","g","r","i","z","y"]
inicios_asoc = {filtro:"filter = '%s'"%(filtro) for filtro in lista_filtros}
inicios_asoc["all"] = ""

def inicio_constraint(filtros_considerados):
    last = " and "            
    if len(filtros_considerados) == 1 :
        if filtros_considerados[0] == "all":
            return ""
        else:
            return inicios_asoc[filtros_considerados[0]]+last
    else:
        lista_inicios_indiv = []
        for x in filtros_considerados:
            lista_inicios_indiv.append(inicios_asoc[x])
            
        return "("+" or ".join(lista_inicios_indiv)+")"+last
        


# ### PARAMETROS USADOS

# In[10]:

def run_por_filtro(filtros_considerados, FBS, pruebas = ""):  #escribir PRUEBAS si se quiere actuar en directorio para tests, si no, dejar en blanco
    print("FBS usado:", FBS)

    metric = NpairsMetric() 
    # ========================= WFD =================================
    constraint1 = inicio_constraint(filtros_considerados)+"note not like '%DD%'"
    wfd_standard = schedUtils.WFD_no_gp_healpixels(64)#, dec_max=2.5, dec_min=-62.5)
    slicer1 = slicers.HealpixSubsetSlicer(64, np.where(wfd_standard==1)[0])#nside = 64, hpid = The subset of healpix id's to use to calculate the metric.
    bundle1 = mb.MetricBundle(metric, slicer1, constraint1)

    # ========================= DDF =================================
    constraint2 = inicio_constraint(filtros_considerados)+"note LIKE '%DD%'"# Esto daba problemas con :"proposalId > 1"
    slicer2 = slicers.HealpixSlicer(nside=64)
    bundle2 = mb.MetricBundle(metric, slicer2, constraint2)

    # # ========================= Other =================================
    # constraint3 = inicio_constraint(filtros_considerados)+"proposalId = 0"
    # slicer3 = slicers.HealpixSlicer(nside=64)
    # bundle3 = mb.MetricBundle(metric, slicer3, constraint3)

    print("==============================================")
    print("constraint WFD:"+constraint1)
    print("constraint DDF:"+constraint2)
    # print("constraint Other:"+constraint3)
    # ### Directorios


    dbDir = './lsst_cadence/FBS_{}/'.format(FBS) # debe correrse en dir contiguo a lsst_cadence
    outDir = '/data/agonzalez/NpairSECONDS/output_FBS_{0}{2}/Npairs_{1}/'.format(FBS, "".join(filtros_considerados), pruebas)
    if not os.path.exists(os.path.abspath(outDir)):
        os.makedirs(os.path.abspath(outDir),exist_ok=True)

    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

    metricDataPath = '/data/agonzalez/NpairSECONDS/output_FBS_{0}{2}/Npairs_{1}/MetricData/'.format(FBS, "".join(filtros_considerados), pruebas)
    if not os.path.exists(os.path.abspath(metricDataPath)):
        os.makedirs(os.path.abspath(metricDataPath),exist_ok=True)

    print("===================================================")
    print("dbDir :",dbDir)
    print("outDir :",outDir)
    print("metricDataPath :",metricDataPath)
    print("===================================================")


    # ### Ahora el bundle group  (SIN considerar noddf)

    dbRuns = show_opsims(dbDir)
    # print(dbRuns[0])
    # In[9]:


    import time

    start = time.perf_counter()

    # dbRuns = ["footprint_big_wfdv1.5_10yrs"]
    cantidad_runs = len(dbRuns)
    n = 1
    for run in dbRuns:
        bDict = {"WFD":bundle1,"DDF":bundle2} #,"Other":bundle3}
        bundle1.setRunName(run)
        bundle2.setRunName(run)
        # bundle3.setRunName(run)
        bgroup = mb.MetricBundleGroup(bDict,                    opSimDbs[run], metricDataPath, resultDbs[run])
        bgroup.runAll()

        ##################### seguimiento tiempo #############
        actual = time.perf_counter()
        tiempo_transcurrido = actual - start
        tiempo_faltante_aprox = tiempo_transcurrido*(cantidad_runs - n)/n
        m,s = divmod(tiempo_transcurrido,60)
        transcurrido_con_formato = "{0:.0f}m:{1:.0f}s".format(m,s)
        m,s = divmod(tiempo_faltante_aprox,60)
        faltante_con_formato = "{0:.0f}m:{1:.0f}s".format(m,s)
        print("________________________________________________________________________")
        print("Terminado dbRun {}/{}, tiempo transcurrido: {}".format(n,cantidad_runs,transcurrido_con_formato))
        print("approx ETA: {}".format(faltante_con_formato))
        print("________________________________________________________________________")
        n +=1

        #####################################################
    print("    =====================================================================================================     ")
    print("    =========================================== FINISHED ================================================     ")
    print("    =====================================================================================================     ")

# In[ ]:

listota_singlelista_filtros = [[x] for x in lista_filtros]
for lista_singlefiltro in listota_singlelista_filtros:
    run_por_filtro(lista_singlefiltro, "1.4.1")




# In[ ]:




