# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 11:31:19 2015

@author: Robert
"""

get_ipython().magic(u'pylab inline')
import scipy
from scipy.integrate import quad, dblquad
import math
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy
from numpy import *
import random

N=100
E0=13
t0=0.5
t=2.
a=math.pow(10,-10)


#Calcul de la densité en dimension supérieure à 3 par la méthode de Monte-Carlo



"""
d=4
def energie_4D(k1,k2,k3,k4):
    return E0-t0-2*t*(cos(k1*a)+cos(k2*a)+cos(k3*a)+cos(k4*a))
    
var=1000000
energie=linspace(E0-t0-2*d*t,E0-t0+2*d*t,100)

cmpt=linspace(0,100,100)
densite=[]
for i in range(100):
    cmpt[i]=0 
 
#En dimension d la zone de Brillouin reste [-pi/a;pi/a]^d

for m in range(var):
    k_1=random.uniform(-pi/a,pi/a)
    k_2=random.uniform(-pi/a,pi/a)
    k_3=random.uniform(-pi/a,pi/a)
    k_4=random.uniform(-pi/a,pi/a)
    E=energie_4D(k_1,k_2,k_3,k_4)
    i=floor((E-(E0-t0-2*d*t))/(4*d*t/100.))
    cmpt[i]+=1
    
 

                    
#On en déduit directement le nombre d'états d'énergie entre E_i et E_(i+1)
for i in range(100):
    densite.append(cmpt[i]/(4*d*t/100.))

"""

"""
d=5
def energie_5D(k1,k2,k3,k4,k5):
    return E0-t0-2*t*(cos(k1*a)+cos(k2*a)+cos(k3*a)+cos(k4*a)+cos(k5*a))
    
var=10000000
energie=linspace(E0-t0-2*d*t,E0-t0+2*d*t,1000)

cmpt=linspace(0,1000,1000)
densite=[]
for i in range(1000):
    cmpt[i]=0 
 
#En dimension d la zone de Brillouin reste [-pi/a;pi/a]^d

for m in range(var):
    k_1=random.uniform(-pi/a,pi/a)
    k_2=random.uniform(-pi/a,pi/a)
    k_3=random.uniform(-pi/a,pi/a)
    k_4=random.uniform(-pi/a,pi/a)
    k_5=random.uniform(-pi/a,pi/a)
    E=energie_5D(k_1,k_2,k_3,k_4,k_5)
    i=floor((E-(E0-t0-2*d*t))/(4*d*t/1000.))
    cmpt[i]+=1
    
 

                    
#On en déduit directement le nombre d'états d'énergie entre E_i et E_(i+1)
for i in range(1000):
    densite.append(cmpt[i]/(4*d*t/1000.))
            
plot(energie,densite)  
"""



d=10
def energie_10D(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10):
    return E0-t0-2*t*(cos(k1*a)+cos(k2*a)+cos(k3*a)+cos(k4*a)+cos(k5*a)+cos(k6*a)+cos(k7*a)+cos(k8*a)+cos(k9*a)+cos(k10*a))
    
var=1000000
energie=linspace(E0-t0-2*d*t,E0-t0+2*d*t,100)

cmpt=linspace(0,100,100)
densite=[]
for i in range(100):
    cmpt[i]=0  
 
#En dimension d la zone de Brillouin reste [-pi/a;pi/a]^d

for m in range(var):
    k_1=random.uniform(-pi/a,pi/a)
    k_2=random.uniform(-pi/a,pi/a)
    k_3=random.uniform(-pi/a,pi/a)
    k_4=random.uniform(-pi/a,pi/a)
    k_5=random.uniform(-pi/a,pi/a)
    k_6=random.uniform(-pi/a,pi/a)
    k_7=random.uniform(-pi/a,pi/a)
    k_8=random.uniform(-pi/a,pi/a)
    k_9=random.uniform(-pi/a,pi/a)
    k_10=random.uniform(-pi/a,pi/a)
    E=energie_10D(k_1,k_2,k_3,k_4,k_5,k_6,k_7,k_8,k_9,k_10)
    i=floor((E-(E0-t0-2*d*t))/(4*d*t/100.))
    cmpt[i]+=1
    


                    
#On en déduit directement le nombre d'états d'énergie entre E_i et E_(i+1)
for i in range(100):
    densite.append(cmpt[i]/(4*d*t/100.))
            
plot(energie,densite)  



