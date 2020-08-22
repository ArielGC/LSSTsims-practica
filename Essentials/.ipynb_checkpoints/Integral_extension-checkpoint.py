import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad, dblquad, cumtrapz, simps
from quasarlf.pubtools.utilities import *
import time

#start = time.time()

def srad(sqdeg): # de sqdeg a srad
	return (np.pi/180)**2 * sqdeg



c = 300000 #km/s
# ========================= Modelo SHELLQS ================================

H0 = 70/10**(-3) #km/s/Gpc      #70/10**6 #km/s/pc
OmegaM = 0.3
OmegaLamb = 0.7
########### ________________ 5.9 < z < 6.5 _________________________

# k = -0.7
# Phi_st = 8.1 # Gpc^-3 mag^-1      #8.1/10**27  # pc^-3 mag^-1
# M_st = -25.3 # Errores : (-1.15, +1.05)
# alpha = -1.39 # Errores : (-0.32,0.45)
# beta = -2.79 # Errores: (-0.48,0.32)

############ ________________ 5.7 < z < 6.5 _________________________
k = -0.7
Phi_st = 10.9 # Gpc^-3 mag^-1      #8.1/10**27  # pc^-3 mag^-1
M_st = -24.9 # Errores : (-0.9, +0.75)
alpha = -1.23 # Errores : (-0.34, 0.44)
beta = -2.73 # Errores: (-0.31, 0.23)
# ==========================================================================
def NtotSimps(my,omega,tupla_z,deltaz1=1/1000,deltaz2=1/1000): #VER ESTOOOOOOOO PORQUE ANTES EL deltaz1 = z/1000 !!!!!!!!!!!!!

	zetta = -2.5*np.log10(1450/9749) #9749 = mean wavelenght of Y filter

	#my = 24 # limiting magnitude in Y band of lsst# esto era para probar en etapas iniciales
	def I(z):
		rango = np.arange(0,z,deltaz1)
		y = 1/np.sqrt(OmegaM*(1+rango)**3+OmegaLamb)
		resultado = simps(y,rango)
		return resultado
		#return quad(lambda zta: 1/np.sqrt(OmegaM*(1+zta)**3+OmegaLamb),0,z)[0] # [0] porque el quad entrega tupla con (resultado,error, otras cosas...)

	def M_lim(z):
		Dabs = 10*10**(-9) #10 #pc
		resultado = my - 5 * np.log10( c*(1+z)*I(z)/(Dabs*H0) ) - (2.5/2) * np.log10(z+1) + zetta
		return resultado


	def dVdz(z):
		return omega * (c**3/H0**3) * (I(z))**2 * 1/np.sqrt(OmegaM*(1+z)**3+OmegaLamb) 

	z_i, z_f = tupla_z 

	rango_z = np.arange(z_i,z_f,deltaz2) #de z =6 a z=14 en vez de inf xd

	""" Para cada z en el rango_z necesito el valor de la integral interior, voy a tener una := lista_exterior. Luego
	se debe usar simps o cumtrapz a esta lista_exterior
	Cada integral interior viene dada por usar simps o cumtrapz en la lista retornada por return_qlf_in_band 
	previamente recortada hasta la M_lim """
	def integral_interior(z):
		M_limite = M_lim(z)
		""" La lista que me entrega return_qlf_in_band probablemente no posee este M_limite en sus valores de M, asi que 
		se interpola este valor phi(M_limite) """
		rangoM,Phi = return_qlf_in_band(z, -5, model='A') #z= 2 , nu = -5 := M1450 # mag, n/mag
		Phi =10**(Phi+9) #Correccion de log10[mag-1Mpc-3] a [mag-1Gpc-3]      
		"""habria que hacer aqui un if si es que se diera milagorsamente el caso que M_lim es justo un valor presente en 
		return_qlf_in_band """
		#print(rangoM)
		#print(M_limite)
		Phi_M_limite = np.interp(M_limite,rangoM,Phi)
		M_menor_aMlim = np.array([x for x in rangoM if x < M_limite])
		Phi_correspondiente = Phi[-len(M_menor_aMlim):]
		# el rangoM que entrega return_qlf_in_band esta al reves e.g. -11, -11.5, -12 ...#
		#print(M_menor_aMlim)        
		M_menorigual_aMlim = np.insert(M_menor_aMlim, 0, M_limite, axis=0)
		Phi_correspondiente = np.insert(Phi_correspondiente,0,Phi_M_limite, axis=0 )
		#print("len(Phi_correspondiente)", len(Phi_correspondiente))
		#print("len(M_menorigual_aMlim)", len(M_menorigual_aMlim))
		return simps(Phi_correspondiente,M_menorigual_aMlim)
# 		return simps(Phi[rangoM<=M_limite],rangoM[rangoM<=M_limite])

#aqui el rango M no es de -200 como antes sino que de -37 masomenos que es casi lo mismo, porque ahi la LF es despreciable 
#(esto a partir mas o menos de M= -30)

	lista_exterior = [integral_interior(z)*dVdz(z) for z in rango_z]
	#end = time.time()
	#delta = end - start
	#print("tiempo (s) = ",delta) #print("took %.2f seconds to process")%(delta)
	#print("len(lista_exterior)", len(lista_exterior))
	return -simps(lista_exterior, rango_z) # el - porque las magnitudes venian al reves entonces al integrar sale un menos

	#return dblquad(lambda M,z: Phi(M,z)*dVdz(z) ,6,np.inf,-200,lambda z :M_lim(z))[0] # -200 porque si pongo -inf  python jode xd


#print("=========================== SHELLQS ==============================")
#print("NQsos  (m_y = 24.4, area = 900 deg2, 5.7 < z < 7) = ",NtotSimps(24.4,srad(900),[5.7,7])) # 0.274 sqrad = 900deg^2
#print ("Descubiertos a 2018 : 93")
#print("")
#end = time.time()
#delta = end - start
#print("tiempo (s) = ",delta) #print("took %.2f seconds to process")%(delta)

#print("=========================== SDSS-Main ==============================")
#print("NQsos  (m_y ~ 20, area = 11240 deg2, 5.7 < z < 6.4))) = ",NtotSimps(20,srad(11240),[5.7,6.4])) # 3.424 sqrad = 11240deg^2
#print ("Descubiertos a 2016 : 29")
#print("")
#print("=========================== SDSS-S82 ==============================")
#print("NQsos  (m_y ~ 22, area = 277 deg2, 5.7 < z < 6.4))) = ",NtotSimps(22,srad(277),[5.7,6.4])) # 3.424 sqrad = 11240deg^2
#print ("Descubiertos a 2016 : 13")
#print("")
#print("=========================== Pan-STARS ==============================")
#print("NQsos  (m_y = 21.1, area = ~2.1*pi rad^2 , 5.6 < z < 6.7))) = ",NtotSimps(21.1,2.1*np.pi,[5.6,6.7])) # 3.424 sqrad = 11240deg^2
#print ("Descubiertos a 2016 : 77")
#print("")
# print("N total de Qsos = ",NtotSimps(24.4,srad(900),[5.7,7]))
#end = time.time()
#delta = end - start
#print("tiempo (s) = ",delta) #print("took %.2f seconds to process")%(delta)
