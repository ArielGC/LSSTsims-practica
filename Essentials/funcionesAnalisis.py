import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pickle

def srad(sqdeg): # de sqdeg a srad
	return (np.pi/180)**2 * sqdeg

def sqdeg(srad):
    return srad*(180/np.pi)**2


def SaveFile(archivo,nombre_guardado):
    a_file = open(nombre_guardado, "wb")
    pickle.dump(diccionario, a_file)
    a_file.close()
    print("Done")

def ResumenMULT_OLD(bDicts,runName,modo,lista_zonas=None):
    if lista_zonas == None:
        lista_zonas = [["ALL","(filter = 'y' or filter = 'z')"],
                       ["WFD","(filter = 'y' or filter = 'z') and note not like 'DD%'"],
                       ["DDF","(filter = 'y' or filter = 'z') and proposalId > 1"]
                      ] ############ OJO CON ESTO HAY QUE CAMBIARLO CUANDO SE ANALICEN SOLO 2 AREASSSSSSSSSSSSSSS !!
    nzonas = len(lista_zonas) # creo que es lo mismo que N_inicial
    bDict = bDicts[runName]
    lista = np.array([x for x in bDict.items()]) 
    # lista es de la forma 
#     [
#     [(1, 'NQsoMetric') <MetricBundle object at 0x7f4842ed54a8>]  ->gatillante1 (i.e. en zona1)
#     [(2, 'NQsoMetric_Area')  <MetricBundle object at 0x7f4842e94ef0>]  -> son1
#     [(3, 'NQsoMetric_y_Ntot')  <MetricBundle object at 0x7f4842dbae10>]  -> son2
#     [(4, 'NQsoMetric_zy_Ntot')  <MetricBundle object at 0x7f4842ed5978>] -> son3
#     .... ]

    a = np.arange(len(lista))
    n = nzonas #N_inicial
    Dgatillantes = 4 #  = distancia de separacion entre metrics gatillantes: [gatillante,son1,son2,son3,gatillante2,son4,son5,son6...etc]
    rango_invalido = np.arange(0,n)*Dgatillantes 
    mask_rango = np.array([False if x in rango_invalido else True for x in a])
    
    array_resultados = lista[mask_rango]
    n_submetrics = Dgatillantes-1
    array_resultados.shape = (int(len(array_resultados)/n_submetrics),n_submetrics,2) 
    # array_resultados tiene esta forma
#     [
#     [ [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] ]  ->zona1
#     [ [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] ]  ->zona2
#     [ [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] ]  ->zona3
#     .... ]   (cada bundle tiene n elem , los que antes eran la key y el value)
    
#     resultado_bundle[0] = [bDictItem1Zx,bDictItem2Zx]      # [0]' porque da lo mismo, ambos de estos bundles tienen igual constraint
#     tal que bDictItem1Zx = [(metricId,metricName),bundle]
#     la shape de resultado_bundle = (len(lista_zonas),2)
    orden = []
    for resultadobundle in array_resultados:
        constraint_este_par = resultadobundle[0][1].constraint
        for zona in lista_zonas:
            constraint_zona = zona[1]
            nombre_zona = zona[0]
            if constraint_este_par == constraint_zona:
                orden.append(nombre_zona)
                
    listaReturn = []
    n=0
    for x in array_resultados: # lista [resulzona1,resulzona2...]
#         print(x[0][0][1])
#         print(x[1][0][1])
#         print(x[2][0][1])
#         print("*******")
        Area = sqdeg(np.sum(x[0][1].metricValues))
        yNtot = np.sum(x[1][1].metricValues)
        zyNtot = np.sum(x[2][1].metricValues)
        nombre_bundle = x[0][0][1]
        ####zona = x[2]
        zona = orden[n]
        #nombre_bundle = x[0][0].split("NQsoMetric10mv2")[1] 
        if modo == "print":
            print("==================================================================================")        
            #print("Proviene de bundle = ",orden[n])
            print("opSim: ",runName)
            print("zona =", zona)
            #print("nombre de este bundle: ",nombre_bundle)
            print("zyNtot = ",round(zyNtot,2))
            print("yNtot = ",round(yNtot,2))
            print("Area tot = ",round(Area,2),"sqdeg")
            print("==================================================================================")
        else:
            listaReturn.append((zyNtot,yNtot,Area,zona,runName)) 
            """
            Ntot = cantidad
            Area = en Sqdeg
            zona = WFD, DDF o ALL
            runName = nombre opsim
            """
        n += 1
    return listaReturn # esta sera tq = [(Ntot1,area1,zona1,runname1),...(totn,arean,zonan,runnamen)]
# _______________________________________________________________________________________________
# _______________________________________________________________________________________________
# ºººººººººººººººººººººººººººººººººººººººººº NPAIRS ººººººººººººººººººººººººººººººººººººººººººººº
# _______________________________________________________________________________________________
# _______________________________________________________________________________________________
def plot_nphist(loghist,View = "xlogy", CutEmptyEdges = False):
    hist, bins = loghist
    print("lenhist",len(hist))
    print("lenbins",len(bins))
    #######################
    if CutEmptyEdges == True:
        nder = 0 # pasos hasta que salga un no nulo (partiendo del primer item) (desde la derecha) i.e indice primer elemento no nulo
        nizq = 0
        for x in hist[::-1]:
            if x == 0:
                nder += 1
            else:
                break
                
        for x in hist:
            if x == 0:
                nizq += 1
            else:
                break
        hist = hist[nizq:len(hist)-nder]
        bins = bins[nizq:len(bins)-nder]
    width = (bins[1:] - bins[:-1])
    center = (bins[:-1] + bins[1:]) / 2
    #####################
    
    plt.bar(center, hist ,align='center', width=width)
    plt.xlabel(r"$\Delta t (dias)$")
    plt.ylabel(r"$N_{pares}$")
    
    if View == "xy":
        return plt.show()
    elif View == "xlogy":
        plt.yscale("log")
        return plt.show()
    elif View == "logxy":
        plt.xscale("log")
        return plt.show()
        
    elif View == "logxlogy":
        plt.xscale("log")
        plt.yscale("log")
        return plt.show()
    
    
def ResumenDEPTH(bDicts,filtro,runName,lista_zonas=None):
    if lista_zonas == None:
        lista_zonas = [["WFD","filter = '{}' and note not like 'DD%'".format(filtro)],
                       ["DDF","filter = '{}' and proposalId > 1".format(filtro)]
                      ]
    nzonas = len(lista_zonas) # creo que es lo mismo que N_inicial
    bDict = bDicts[runName]
    lista = np.array([x for x in bDict.items()]) 
    # lista es de la forma 
#     [
#     [(1, 'Nombre de metric'), <MetricBundle object at 0x7f4842ed54a8>] ->zona1
#     [(2, 'Nombre de metric'), <MetricBundle object at 0x7f4842ed54a8>] ->zona2 
#               ]  
# En este no deberia ser necesario verificar orden porque no se crean nuevos bundles, solo se
# actualizan los ya creados, entonces es el mismo orden que yo le di, i.e. WFD, DDF. Pero igual , 
# por siaca
    orden = []
    for resultadobundle in lista:
        constraint_este_par = resultadobundle[1].constraint
        for zona in lista_zonas:
            constraint_zona = zona[1]
            nombre_zona = zona[0]
            if constraint_este_par == constraint_zona:
                orden.append(nombre_zona)
                
    listaReturn = []
    n=0
    for x in lista: 
        coadded_max = np.max(x[1].metricValues)
        coadded_min = np.min(x[1].metricValues)
        nombre_bundle = x[0][1]
        zona = orden[n]
        
        listaReturn.append(
            (coadded_min,coadded_max,zona,runName)  )
        n += 1
    return listaReturn # esta sera tq = [(coadded_minn,coadded_maxn,zonan,runNamen),...(coadded_minn,coadded_maxn,zonan,runNamen)]

class HiddenPrints: # pa no printear si se pone with <esto>: <codigo>
    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, 'w')

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout
        
def ResumenMULT_NEW(bDict,runName,filtros,modo,lista_zonas=None):
    """
    bDict.items() es de la forma 
    ((1, 'NtotMetric'), <Bundle>)---------------------   >>>>>>>> gatillante 1
    ((2, 'NtotMetric_Area'), <Bundle>)               |
    ((3, 'NtotMetric_f1f2o10_Ntot'), <Bundle>)       |
    ((4, 'NtotMetric_f1f2o5_Ntot'), <Bundle>)        |------------ ZONA 1
    ((5, 'NtotMetric_f2o10_Ntot'), <Bundle>)         |
    ((6, 'NtotMetric_f2o5_Ntot'), <Bundle>)-----------
    ((7, 'NtotMetric'), <Bundle>)---------------------   >>>>>>>> gatillante 2
    .                                                |
    .                                                |------------ ZONA 2
    .                                                |
    ((12, 'NtotMetric_f2o5_Ntot'),<Bundle>)-----------
    """
    f1 = filtros[0]
    f2 = filtros[1]
    if lista_zonas == None:
        lista_zonas = [
                        ["WFD", "(filter = '{}' or filter = '{}') and note not like 'DD%'".format(f1, f2)],
                        ["DDF", "(filter = '{}' or filter = '{}') and proposalId > 1".format(f1, f2)]
                      ] 
    
    nzonas = len(lista_zonas) # creo que es lo mismo que N_inicial
    lista = np.array([x for x in bDict.items()]) 

    a = np.arange(len(lista))
    n = nzonas #N_inicial
    Dgatillantes = 6 #4 #  = distancia de separacion entre metrics gatillantes (restar sus indices): [gatillante,son1,son2,son3,gatillante2,son4,son5,son6...etc] Dgat = 4 
    rango_invalido = np.arange(0,n)*Dgatillantes 
    mask_rango = np.array([False if x in rango_invalido else True for x in a])
    array_resultados = lista[mask_rango]   # -> se eliminan los gatillantes de la lista
    
    n_submetrics = Dgatillantes-1
    array_resultados.shape = (int(len(array_resultados)/n_submetrics),n_submetrics,2) 
    lista_resultados_por_zona = array_resultados
    # lista_resultados_por_zona tiene ahora esta forma
#     [
#        [     [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] etc...    ]  ->zona1
#        [     [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] etc...    ]  ->zona2
#        [     [(id,key),bundleArea], [(id,key),Bundle yNtot], [(id,key),bundle zyNtot] etc...    ]  ->zona3
#  ...]   (cada bundle tiene n elem , los que antes eran la key y el value)
    
#     resultado_bundle[0] = [bDictItem1Zx,bDictItem2Zx]      # [0]' porque da lo mismo, ambos de estos bundles tienen igual constraint
#     tal que bDictItem1Zx = [(metricId,metricName),bundle]
#     la shape de resultado_bundle = (len(lista_zonas),2)
    orden = []
    for resultados_zona in lista_resultados_por_zona:
        constraint_esta_zona = resultados_zona[0][1].constraint
        for zona in lista_zonas:
            constraint_zona = zona[1]
            nombre_zona = zona[0]
            if constraint_esta_zona == constraint_zona:
                orden.append(nombre_zona)
    print("orden:", orden)
                
    listaReturn = []
    n=0
    for x in lista_resultados_por_zona: # lista [resulzona1,resulzona2...]
#         print(x[0][0][1])
#         print(x[1][0][1])
#         print(x[2][0][1])
#         print("*******")
        Area =   sqdeg(np.sum(x[0][1].metricValues))
        Ntot_f1f2o10 = np.sum(x[1][1].metricValues)
        Ntot_f1f2o5 =  np.sum(x[2][1].metricValues)
        Ntot_f2o10 =   np.sum(x[3][1].metricValues)
        Ntot_f2o5 =    np.sum(x[4][1].metricValues)

        nombre_bundle = x[0][0][1]
        ####zona = x[2]
        zona = orden[n]
        #nombre_bundle = x[0][0].split("NQsoMetric10mv2")[1] 
        if modo == "print":
            print("==================================================================================")        
            #print("Proviene de bundle = ",orden[n])
            print("opSim: ",runName)
            print("zona =", zona)
            #print("nombre de este bundle: ",nombre_bundle)
            print("Ntot_{}{}o10 =".format(f1,f2),Ntot_f1f2o10)
            print("Ntot_{}{}o5 =".format(f1,f2),Ntot_f1f2o5)
            print("Ntot_{}o10 =".format(f2),Ntot_f2o10)
            print("Ntot_{}o5 =".format(f2),Ntot_f2o5)  
            print("Area tot = ",round(Area,2),"sqdeg")
            print("==================================================================================")
        else:
            listaReturn.append((Ntot_f1f2o10,Ntot_f1f2o5,Ntot_f2o10,Ntot_f2o5,Area,zona,filtros,runName)) 
            """
            Ntot = cantidad
            Area = en Sqdeg
            zona = WFD, DDF o ALL
            runName = nombre opsim
            """
        n += 1
    return listaReturn # esta sera tq = [(Ntot1,area1,zona1,runname1),...

def Histogramita(lista_listasreturn,zona_interes,variante):    
    #lista_Ntot = []
    lista_NtArRn = [] #Ntot, area y Runname
    for lista_return in lista_listasreturn:
        print(lista_return)
        for tupla in lista_return: # en cada tupla resultado se selecciona unicamente el de la zona deseada
            Ntot_f1f2o10,Ntot_f1f2o5,Ntot_f2o10,Ntot_f2o5,Area,zona,filtros,runName = tupla
            DictNtot = {"f1f2o10":Ntot_f1f2o10,"f1f2o5":Ntot_f1f2o5,"f2o10":Ntot_f2o10,"f2o5":Ntot_f2o5}
            if zona == zona_interes:
                lista_NtArRn.append((round(DictNtot[variante],2),round(Area,2),runName)) 
                break
    
    # filtros quedo definido como en el ultimo for, pero todos tienen los mismo filtros asi que da igual    
    f1,f2 = filtros
    x = np.arange(len(lista_NtArRn)) # ni lo ocupo xd
    lista_NtArRn = sorted(lista_NtArRn,key=lambda x:x[0])
    y = [i[0] for i in lista_NtArRn]


    plt.hist(y,color ="gray")
    plt.title(zona_interes)
    plt.ylabel("Cantidad Opsims")
    plt.xlabel("Ntot_{}".format(variante))
    plt.show()
    print("formato :      (Ntot_{}{}o{}   ,Area(sqdeg)    ,opsim)".format(f1,f2,variante.split("o")[1]))
    print("Max Ntot_{}{}o{} = ".format(f1,f2,"o"+variante.split("o")[1]), max(lista_NtArRn,key=lambda x:x[0]))
    print("2nd Max Ntot_{}{}o{} = ".format(f1,f2,variante.split("o")[1]), lista_NtArRn[-2])
    print("Min Ntot_{}{}o{} = ".format(f1,f2,variante.split("o")[1]),min(lista_NtArRn,key=lambda x:x[0]))
    print("Mediana Ntot_{}{}o{} =".format(f1,f2,variante.split("o")[1]), round(np.median(y),2))
    print("Media Ntot_{}{}o{} = ".format(f1,f2,variante.split("o")[1]), round(np.mean(y),2))
    print("StD = ",round(np.std(y),2))
    

def Resumen_NtotV2(runNames,bundleDicts):    
    lista_dict_simplificados = [] #opsim,zona,variante1:values1,variante2:values2... -> opsim:footprint2_psa,zona:"WFD","g_o10_Ntot":val...
    for runName in runNames:
        bDict = bundleDicts[runName]
        ## bDict simplificado: {""}
        bDict_simplificado_WFD = {"opsim":runName, "zona":"WFD"}
        bDict_simplificado_DDF = {"opsim":runName, "zona":"DDF"}

        for item in bDict.items(): #   item de la forma ((3, 'NtotMetricV2_g_o10_Ntot'), <lsst.sims.maf.metricBundles.metricBundle.MetricBundle object at 0x7fee4a609358>)

            if item[0][1] == "NtotMetricV2":
                continue
            else:
                variante = item[0][1][13:] # elimina el NtotMetricV2
                if variante[-4:].lower() == "area": variante = "Area"
                values = round(np.sum(item[1].metricValues),2)

            if item[1].constraint == "note not like 'DD%'": # item[1] es el  <bundle>
                bDict_simplificado_WFD[variante] = values
            elif item[1].constraint == "proposalId > 1":
                bDict_simplificado_DDF[variante] = values
            else: raise ValueError("no pertenece a ninguna zona wtf")
        if len(bDict_simplificado_DDF) > 2: lista_dict_simplificados.append(bDict_simplificado_DDF)
        if len(bDict_simplificado_WFD) > 2: lista_dict_simplificados.append(bDict_simplificado_WFD)
    return lista_dict_simplificados

def Resumen_Ntot_intV5(runNames,bundleDicts):     # PARA SER USADO CON METRICS CORRIDAS CON V5 LOOKUPTABLES, evitar prob con bigwfd (%DD%)
    lista_dict_simplificados = [] #opsim,zona,variante1:values1,variante2:values2... -> opsim:footprint2_psa,zona:"WFD","g_o10_Ntot":val...
    for runName in runNames:
        bDict = bundleDicts[runName]
        ## bDict simplificado: {""}
        bDict_simplificado_WFD = {"opsim":runName, "zona":"WFD"}
        bDict_simplificado_DDF = {"opsim":runName, "zona":"DDF"}

        for item in bDict.items(): #   item de la forma ((3, 'NtotMetricV2_g_o10_Ntot'), <lsst.sims.maf.metricBundles.metricBundle.MetricBundle object at 0x7fee4a609358>)

            if item[0][1] == "NtotMetricV2":
                continue
            else:
                variante = item[0][1][13:] # elimina el NtotMetricV2
                if variante[-4:].lower() == "area": variante = "Area"
                values = round(np.sum(item[1].metricValues),2)
                
            if item[1].constraint == "note NOT LIKE '%DD%'": # item[1] es el  <bundle>
                bDict_simplificado_WFD[variante] = values
            elif item[1].constraint == "note LIKE '%DD%'":
                bDict_simplificado_DDF[variante] = values
            else: raise ValueError("no pertenece a ninguna zona wtf")
            
        if len(bDict_simplificado_DDF) > 2: lista_dict_simplificados.append(bDict_simplificado_DDF)
        if len(bDict_simplificado_WFD) > 2: lista_dict_simplificados.append(bDict_simplificado_WFD)
    return lista_dict_simplificados
