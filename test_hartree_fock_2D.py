
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:23:55 2015

@author: Robert
"""
import numpy
from numpy import *
import scipy
from scipy.integrate import quad, dblquad
import math
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D

e2=2.3*math.pow(10,-28)

#Integrale en jeu dans le terme de Hartree Fock en 2D

N=2
a=math.pow(10,-10)
d=a #etalement caractéristique de la fonction d'onde orbitale

k_x=linspace(-pi/a,pi/a,N)
k_y=linspace(-pi/a,pi/a,N)

#Gaussienne
def chi(x,y):
    return exp(-(x**2+y**2)/(d/a)**2)
  
def integrande_1(l1,p1,s,t,u,v):
    return (chi(u-l1,v-p1)**2)/sqrt((u-s)**2+(v-t)**2)

def integrande_2(l1,p1,s,t):
   return dblquad(lambda u, v : integrande_1(l1,p1,s,t,u,v), p1-5*(d/a), t, lambda u : l1-5*(d/a), lambda u : s )[0]+dblquad(lambda u, v : integrande_1(l1,p1,s,t,u,v), t, p1+5*(d/a), lambda u : l1-5*(d/a), lambda u : s )[0]+dblquad(lambda u, v : integrande_1(l1,p1,s,t,u,v), p1-5*(d/a), t, lambda u : s, lambda u : l1+5*(d/a) )[0]+dblquad(lambda u, v : integrande_1(l1,p1,s,t,u,v), t, p1+5*(d/a), lambda u : s, lambda u : l1+5*(d/a) )[0]


#Facteur a**3 près
def I_2D(l1,p1,l2,p2):
    return dblquad(lambda s, t : (chi(s-l2,t-p2)**2)*integrande_2(l1,p1,s,t),p2-5*(d/a), p2+5*(d/a), lambda s : l2-5*(d/a), lambda s : l2+5*(d/a))[0]
   
   
print(I_2D(0.,0.,0.,0.))