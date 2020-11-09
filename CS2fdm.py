#-*- coding: utf-8 -*-
"""
Created on Wed Feb 12 14:35:01 2020

author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
from FDM1D import *
import matplotlib.pyplot as plt

l=1
numelem=50
ite=10

#Problem DATA
#w=2*3.141592/period
w=10000

#frequency [Hz]
#lamb=(w*8.08*(10**-3))**(-0.5)
#BC
T0=1
dq=l/numelem
s0=0.5

#Setting data
L=initL(numelem,l)
T=initT(numelem,T0)

for i in range(ite):
    #Mechanical part
    K=initK(numelem,dq,w,T)
    F=initF(numelem,s0)

    S=np.linalg.solve(K,F)
    
    #Getting Stresses
    S1=initS1(numelem,s0)
    S2=initS2(numelem)
    S1[1:-1]=S[:numelem-1]
    S2[1:-1]=S[numelem-1:]

    #Thermal part
    FT=initFT(numelem,dq,S1,S2,T)
    KT=initKT(numelem,dq)
    T[:-1]=np.linalg.solve(KT,FT)


plt.plot(L,S1)
plt.plot(L,S2)
#T=T*190-125
#T=5*(T-32)/9
#plt.plot(L,T)
plt.grid()


mukx1 = [0,0.0243,0.1096,0.1802,0.2290,0.2777,0.3290,0.3974,0.4585,0.5124,0.5699,0.6250,0.6789,0.7352,0.7928,0.8504,0.9226,0.9874] 
mukS1 = [0,-0.0400,-0.1557,-0.2415,-0.2972,-0.3333,-0.3502,-0.3289,-0.2829,-0.2268,-0.1418,-0.0566,0.0383,0.1283,0.2377,0.3276,0.4169,0.4676]

mukx2 = [0,0.0293,0.0635,0.2356,0.3210,0.4052,0.5101,0.5809,0.6565,0.7628,0.8483,0.9058,0.95349,0.9987]
mukS2 = [0,0.0132,0.0117,-0.0249,-0.0530,-0.0859,-0.1294,-0.1617,-0.1796,-0.1794,-0.1394,-0.0981,-0.0564,-0.0048]
plt.scatter(mukx1,mukS1)
plt.scatter(mukx2,mukS2)
