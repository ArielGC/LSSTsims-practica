import numpy as np 
import matplotlib.pyplot as plt 
from scipy.integrate import quad, dblquad
import decimal


c = 300000 #km/s
## _____________ Modelo SHELLQS _________
H0 = 70/10**(-3) #km/s/Gpc      #70/10**6 #km/s/pc
OmegaM = 0.3
OmegaLamb = 0.7
k = -0.7
Phi_st = 8.1 # Gpc^-3 mag^-1      #8.1/10**27  # pc^-3 mag^-1
M_st = -25.3 # Errores : (-1.15, +1.05)
alpha = -1.39 # Errores : (-0.32,0.45)
beta = -2.79 # Errores: (-0.48,0.32)
def Ntot(my,omega):
	def Phi(M,z):
		denominador = 10**(   0.4*(alpha + 1)*(M - M_st)   ) + 10**(   0.4*(beta + 1)*(M - M_st)   )
		return 10**(k*(z-6)) * Phi_st / denominador
	##________ Porcion cielo : _______________
	#omega = 4*np.pi# esto era para probar en etapas iniciales
	##________________________________________ 
	zetta = -2.5*np.log10(1450/9749) #9749 = mean wavelenght of Y filter

	#my = 24 # limiting magnitude in Y band of lsst# esto era para probar en etapas iniciales
	def I(z):
		return quad(lambda zta: 1/np.sqrt(OmegaM*(1+zta)**3+OmegaLamb),0,z)[0] # [0] porque el quad entrega tupla con (resultado,error, otras cosas...)

	def M_lim(z):
		Dabs = 10*10**(-9) #10 #pc
		return my - 5 * np.log10( c*(1+z)*I(z)/(Dabs*H0) ) - (2.5/2) * np.log10(z+1) + zetta


	def dVdz(z):
		return omega * (c**3/H0**3) * (I(z))**2 * 1/np.sqrt(OmegaM*(1+z)**3+OmegaLamb) 

	return dblquad(lambda M,z: Phi(M,z)*dVdz(z) ,6,np.inf,-200,lambda z :M_lim(z))[0] # -200 porque si pongo -inf  python jode xd
