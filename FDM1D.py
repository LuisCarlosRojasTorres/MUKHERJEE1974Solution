# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 11:15:40 2020

author: Luis Carlos A. Rojas Torres
email: luiscarlos.bsf@oceanica.ufrj.br
"""
import numpy as np
import matplotlib.pyplot as plt

def initL(numelem,l):
    L=np.linspace(0,l,numelem+1)
    return L

def initT(numelem,T0):
    T=np.zeros(numelem+1)
    for i in range(T.size):
        T[i]=T0
    return T
   
def d1(dq,w,T):
    # calcula dq**2 * d1 
    dq2=dq**2
    beta=-0.214
    c1=4.61*(10**-11)
    gamma=3.21
    #a = l
    a=1.023*(10**-3)*(w**(2+beta))*190**gamma
    return dq2*c1*a*T**gamma
    

def d2(dq,w,T):
    # calcula dq**2 * d1
    dq2=dq**2
    beta=-0.214
    c2=1.62*(10**-11)
    gamma=3.21
    a=1.023*(10**-3)*(w**(2+beta))*190**gamma
    return dq2*c2*a*T**gamma
    

def initK(numelem,dq,w,T):
    K=np.zeros((2*(numelem-1),2*(numelem-1)))
    for i in range(numelem-1):
        K[i][i]=d1(dq,w,T[i+1])-2
        K[i][numelem-1+i]=d2(dq,w,T[i+1])
        
        K[numelem-1+i][numelem-1+i]=d1(dq,w,T[i+1])-2
        K[numelem-1+i][i]=-d2(dq,w,T[i+1])
    
    for i in range(numelem-2):
        K[i][i+1]=1
        K[i+1][i]=1
        K[numelem-1+i][numelem+i]=1
        K[numelem+i][numelem-1+i]=1
    return K

def initS1(numelem,s0):
    s1=np.zeros(numelem+1)
    s1[-1]=s0
    return s1

def initS2(numelem):
    s2=np.zeros(numelem+1)
    return s2

def initF(numelem,s0):
    F=np.zeros(2*(numelem-1))
    F[numelem-2]=-s0
    return F

#THERMAL
def initFT(numelem,dq,S1,S2,T):
    FT=np.zeros(numelem)
    for i in range(numelem):
        FT[i]=S2[i]*(S1[i+1]-S1[i])-S1[i]*(S2[i+1]-S2[i])-dq
    FT[numelem-1]-=1
    return FT
    
def initKT(numelem,dq):
    KT=np.zeros((numelem,numelem))
    for i in range(numelem-1):
        KT[i][i]=-1
        KT[i][i+1]=1
        KT[i][0]-=dq
    KT[numelem-1][numelem-1]=-1
    KT[numelem-1][0]=-dq
    return KT    
