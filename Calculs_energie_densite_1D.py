# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:33:59 2015

@author: Robert
"""


get_ipython().magic(u'pylab inline')

import numpy
import scipy
from scipy.integrate import quad, dblquad
import math
from math import *

N=100
E0=13
t0=0.5
t=2
a=math.pow(10,-10)

# Energie pour le réseau d'atomes en 1D avec un pas de a
energie_1D=[(E0-t0-2*t*cos((2*pi*i/(N*a))*a)) for i in range(-N/2,N/2,1)]
abscisses=[(2*pi*i/(N*a)) for i in range(-N/2,N/2,1)]

plot(abscisses,energie_1D)
#legend(['E(k)'])


# In[8]:

E = linspace(E0-t0-2*t,E0-t0+2*t,100)#Dernier argument : nombre de points découpant l'intervalle.

#Densité d'états en 1D

#z = N/(2*t)*(1/sqrt((1-((E-(E0-t0))/(2*t))**2)))

def d_1D(x):
    return N/(2*pi*t)*(1/sqrt((1-((x-(E0-t0))/(2*t))**2)))

z=[]
for i in range(1,99):
    z.append(d_1D(E[i]))

#plt.plot(E[1:99], z,label= "D(E)")
#xlabel("Energie E (eV)")
#ylabel("D(E)")


#Normalisation :
print(quad(lambda x : d_1D(x),E0-t0-2*t,E0-t0+2*t)[0])
#OK : on trouve bien N=100, le nombre total d'états

#print(numpy.arccos((E-(E0-t0))/(2*t)))

nb_etats_1D=N-(N/pi)*numpy.arccos((E-(E0-t0))/(2*t)) # Nombre d'états d'énergie inférieure ou égale à E
#print(nb_etats_1D)

plot(E,nb_etats_1D)
xlabel("Energie E (eV)")
ylabel("N_1D(E)")
