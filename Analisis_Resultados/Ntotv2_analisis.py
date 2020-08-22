import glob
import os
import sys
import pickle
#print (os.getcwd())
sys.path.append("..")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# For the WFD footprint definitions
from lsst.sims.featureScheduler import utils as schedUtils
#print(schedUtils.__file__)
# To calculate metrics with MAF
import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as mb
#print(mb.__file__)
import lsst.sims.maf.plots as plots
import lsst.sims.maf.batches as batches
# import convenience functions

from Essentials.opsimUtils import *
from Essentials.funcionesAnalisis import Resumen_Ntot_intV5, HiddenPrints

# FBS = "1.5"
# mod = "A"
for FBS in ["1.4.1","1.5"]:
    for mod in ["A","B"]:
        # user provided paths
        resultDbPath =  '/data/agonzalez/output_FBS_{}/bolNtot_mod{}_FINAL/'.format(FBS,mod)
        metricDataPath = '/data/agonzalez/output_FBS_{}/bolNtot_mod{}_FINAL/MetricData/'.format(FBS,mod)
        # get a dictionary of resultDb from given directory
        resultDbs = getResultsDbs(resultDbPath)
        #print(resultDbs)
        # the following line will be useful if you did not run MAF on all 75 opsims
        runNames = list(resultDbs.keys())

        # for key in resultDbs.keys():
        #     print(resultbDbs[key].getMetricDisplayInfo()) #for x in resultDbs: print(x)
        bundleDicts = {}
        with HiddenPrints():
            for runName in resultDbs:
                bundleDicts[runName] = bundleDictFromDisk(resultDbs[runName], runName, metricDataPath)

###################################### ACCION #########################################
        lista_dict_simplificados = Resumen_Ntot_intV5(runNames,bundleDicts)

        archivo = open("lista_dict_simplificados{}{}.pkl".format(FBS,mod),"wb")
        pickle.dump(lista_dict_simplificados,archivo)
        archivo.close()







