#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import inspect
import glob
import sys
# sys.path.append("/home/idies/workspace/Storage/aagonzalez6/persistent/LSST_OpSim/Scripts_NBs/Essentials")

import numpy as np
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
#%matplotlib inline #con nbconvert esto daba problema

from Essentials.FMNtot import NtotMetricV2

inspect.getfile(NtotMetricV2)


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


# In[ ]:





# ### PARAMETROS USADOS

# In[3]:


def ResultadosNtotBolV2(FBS, mod):

    # ==========================================================
#     mod = "A"
    # FBS = "1.5"
    # modo = "A"
    # filtros_considerados = ["u","g"]  # f1,f2 tq f2 mas rojo que f1 
    # ==========================================================
    #validacion(filtros_considerados)
    #f1,f2 = filtros_considerados
    # g_modA_LookupT_extension.pk
#     lookup_table = "{}_mod{}_LookupT_extension.pkl".format(f2, modo) # debe estar en la carpeta de /lookuptables en /essentials
                                                                     # f2 porque ese se ocupa , el f1 es para potencial lyman pbreak nomas
    #filtros_modo = "{}_mod{}".format("".join(filtros_considerados),modo)
    print("FBS usado:", FBS)
    print("mod:",mod)

    #####################################################################################
    ################################## 3 BUNDLES ########################################
    #####################################################################################
    metric = NtotMetricV2(mod, f1f2diff = 2) 
    # ========================= WFD =================================
    constraint1 = "note NOT LIKE '%DD%'"
    wfd_standard = schedUtils.WFD_no_gp_healpixels(64)#, dec_max=2.5, dec_min=-62.5)
    slicer1 = slicers.HealpixSubsetSlicer(64, np.where(wfd_standard==1)[0])#nside = 64, hpid = The subset of healpix id's to use to calculate the metric.
    bundle1 = mb.MetricBundle(metric, slicer1, constraint1)

    # ========================= DDF =================================
    constraint2 = "note LIKE '%DD%'"
    slicer2 = slicers.HealpixSlicer(nside=64)
    bundle2 = mb.MetricBundle(metric, slicer2, constraint2)

    print("==============================================")
    print("constraint WFD:"+constraint1)
    print("constraint DDF:"+constraint2)

    #####################################################################################
    ################################# DIRECTORIOS #######################################
    #####################################################################################

    #Please enter your SciServer username between the single quotes below!
   # your_username = 'aagonzalez6'
    # Check avaiable database directoies
    show_fbs_dirs()
   # if your_username == '': # do NOT put your username here, put it in the cell at the top of the notebook.
   #     raise Exception('Please provide your username!  See the top of the notebook.')

    dbDir = './lsst_cadence/FBS_{}/'.format(FBS)
    outDir = '/data/agonzalez/output_FBS_{}/bolNtot_mod{}_FINAL/'.format(FBS,mod)
    if not os.path.exists(os.path.abspath(outDir)):
        os.makedirs(os.path.abspath(outDir),exist_ok=True)

    opSimDbs, resultDbs = connect_dbs(dbDir, outDir)

    metricDataPath = '/data/agonzalez/output_FBS_{}/bolNtot_mod{}_FINAL/MetricData/'.format(FBS,mod)
    if not os.path.exists(os.path.abspath(metricDataPath)):
        os.makedirs(os.path.abspath(metricDataPath),exist_ok=True)

    print("===================================================")
    print("dbDir :",dbDir)
    print("outDir :",outDir)
    print("metricDataPath :",metricDataPath)
    print("===================================================")

    #####################################################################################
    ################################# BUNDLE GROUP ######################################
    #####################################################################################

    dbRuns = show_opsims(dbDir)
    print(dbRuns)
    dbRuns = [x for x in dbRuns if "noddf" not in x]
    
    #archivo70plus = open("jhu70plus{}.txt".format(FBS),"r")
    #dbRuns = [x.rstrip() for x in list(archivo70plus)]
    #archivo70plus.close()
       
    for run in dbRuns: #[70:]:
        bDict = {"WFD":bundle1,"DDF":bundle2}
        bundle1.setRunName(run)
        bundle2.setRunName(run)
        bgroup = mb.MetricBundleGroup(bDict, opSimDbs[run], metricDataPath, resultDbs[run])
        bgroup.runAll()
        #print(bDict)
        #printResumen(bDict)


# In[4]:


import time
class cronometro():
	def __init__(self):
		self.times = [0,0]

	def clicki(self):
		print("__________________________________")
		self.times[0] = time.perf_counter()
	def clickf(self):
		self.times[1] = time.perf_counter()
		penultimo_click = self.times[-2]
		ultimo_click = self.times[-1]
		print("dt 2 ultimos clicks:",ultimo_click-penultimo_click)
		print("     ^^^^^^^^^^^^^^^^^^^^^^^     \n")

cron = cronometro()

listilla_FBS = ["1.5"]
listilla_modos = ["B"]
listilla_parejas_filtros = [["u","g"],["g","r"],["r","i"],["i","z"],["z","y"]]

cron.clicki()
for FBS in listilla_FBS:
    for mod in listilla_modos:
        ResultadosNtotBolV2(FBS, mod)
cron.clickf()


# In[ ]:





# In[ ]:




