# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 18:10:19 2015

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
 

#Calcul de la densité des états d'énergie pour le graphène par la méthode de Monte-Carlo



def energie_graphene_liant(kx,ky):
    return E0-t0-t*sqrt(3 + 2*(2*cos((3*a/2)*kx)*cos((sqrt(3)*a/2)*ky)+cos(sqrt(3)*a*ky)))
def energie_graphene_antiliant(kx,ky):
    return E0-t0+t*sqrt(3 + 2*(2*cos((3*a/2)*kx)*cos((sqrt(3)*a/2)*ky)+cos(sqrt(3)*a*ky)))
    
var=1000000
energie=linspace(E0-t0-3*t,E0-t0+3*t,500)

cmpt=linspace(0,500,500)
densite=[]
for i in range(500):
    cmpt[i]=0 
 
#La zone de Brillouin du graphène est un hexagone inclus dans [-pi/a;pi/a]^2
#Il faut faire des tirages aléatoires de (kx,ky) uniquement dans la 1 Z.B. pour obtenir le bon résultat de densité
#Sinon des points k différents et pas dans la même cellule unité pourraient être tirer et contribuer à la même énergie


#On travaille dans la base (a1*,a2*) du réseau de Bravais réciproque :
A_1=[2*pi/3,-2*pi/sqrt(3)]
A_2=[2*pi/3,2*pi/sqrt(3)]

for m in range(var):
    alpha=random.uniform(0.,1.)
    beta=random.uniform(0.,1.)
    #k1=[(alpha*A_1[i]+beta*A_2[i]) for i in range(2)]
    
    #Il y a deux états d'énergie pour le couple (k_x,k_y)
    #E1=energie_graphene_liant(k1[0],k1[1])
    E1=E0-t0-t*sqrt(3+2*(cos(2*pi*alpha)+cos(2*pi*beta)+cos(2*pi*(beta-alpha))))
    #print(E1)
    ind1=floor((E1-(E0-t0-3*t))/(6*t/500.))
    cmpt[ind1]+=1
    #E2=energie_graphene_antiliant(k1[0],k1[1])
    E2=E0-t0+t*sqrt(3+2*(cos(2*pi*alpha)+cos(2*pi*beta)+cos(2*pi*(beta-alpha))))
    #print(E2)
    ind2=floor((E2-(E0-t0-3*t))/(6*t/500.))
    cmpt[ind2]+=1
    
                    
#On en déduit directement le nombre d'états d'énergie entre E_i et E_(i+1)
for i in range(500):
    densite.append(cmpt[i]/(6*t/500.))


plot(energie,densite)