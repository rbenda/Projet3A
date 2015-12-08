# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:43:17 2015

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
from Calculs_energie_densite_2D import D_2D

N=100
E0=13
t0=0.5
t=2.
a=math.pow(10,-10)

def alpha(E):
    return (E0-t0-E)/(2*t)

#Nombre d'états en 3D
"""
#Il faut appeler la fonction N2D de l'autre code
def N_3D(E):
    if E > E0-t0+2*t and E<=E0-t0+6*t:
        return (N**3/pi)*numpy.arccos(alpha(E)+2)+(N/pi)*quad(lambda z: N_2D(E+2*t*cos(z)),numpy.arccos(alpha(E)+2),pi)[0]
    if E < E0-t0-2*t and E>E0-t0-6*t:
        return (N/pi)*quad(lambda z: N_2D(E+2*t*cos(z)),0,numpy.arccos(alpha(E)-2))[0]
    if E > E0-t0-2*t and E < E0-t0+2*t:
        return (N/pi)*(quad(lambda z: N_2D(E+2*t*cos(z)),0,numpy.arccos(alpha(E)))[0]+quad(lambda z: N_2D(E+2*t*cos(z)),numpy.arccos(alpha(E)),pi)[0])
    if (E<E0-t0-6*t):
        return 0
    if (E>E0-t0+6*t):
        return N_3D(E0-t0+6*t)



nb2=50.


E_3D = np.linspace(E0-t0-7*t,E0-t0+7*t,50)
nb_etats_3D= [N_3D(E0-t0-7*t+i*(14*t)/nb2) for i in range(0,50,1)]
"""
#Visualisation continuité de la pente en E0-t0-6t

#E_3D = linspace(E0-t0-6.5*t,E0-t0-5.5*t,50)
#nb_etats_3D= [N_3D(E0-t0-6.5*t+i*(1*t)/nb2) for i in range(0,200,1)]


#On peut renormaliser par N**3 pour comparer l'évolution du nombre d'états en 2D et en 3D ?
#plot(E_3D,nb_etats_3D)

#plot(E_2D,nb_etats_2D)


# In[39]:

#Densité d'états en 3D
def D_3D(E):
    if (E > E0-t0+2*t):
        return (N/pi)*quad(lambda z: D_2D(E+2*t*cos(z)),numpy.arccos(alpha(E)+2),pi)[0]
    if (E < E0-t0-2*t):
        return (N/pi)*quad(lambda z: D_2D(E+2*t*cos(z)),0,numpy.arccos(alpha(E)-2))[0]
    if E > E0-t0-2*t and E < E0-t0+2*t:
        return (N/pi)*(quad(lambda z: D_2D(E+2*t*cos(z)),0,numpy.arccos(alpha(E)))[0]+quad(lambda z: D_2D(E+2*t*cos(z)),numpy.arccos(alpha(E)),pi)[0])+(N/(2*pi*t))*N_2D(E0-t0-0.000001)/(sqrt(1-alpha(E)**2))
        #return (N/pi)*(quad(lambda z: D_2D(E+2*t*cos(z)),0,numpy.arccos(alpha(E))-0.01*t)[0]+quad(lambda z: D_2D(E+2*t*cos(z)),numpy.arccos(alpha(E))+0.01*t,pi)[0])+(N/(2*pi*t))*N_2D(E0-t0-0.000001)/(sqrt(1-alpha(E)**2))
        #return (N/pi)*(quad(lambda z: D_2D(E+2*t*cos(z)),0,pi))[0]


#densite_etats_3D=[(N_3D(E0-t0-6*t+(i+1)*(12*t)/nb2)-N_3D(E0-t0-6*t+i*(12*t)/nb2))/(12*t/nb2) for i in range(0,nb2-1,1)]


E_3D = linspace(E0-t0-6*t,E0-t0+6*t,50)
densite_etats_3D= [D_3D(E0-t0-6*t+i*(12*t)/50.) for i in range(0,50,1) ]

plt.plot(E_3D,densite_etats_3D)
"""
#Calcul de la densité en 3D par la méthode des histogrammes

def energie_3D(kx,ky,kz):
    return E0-t0-2*t*(cos(kx*a)+cos(ky*a)+cos(kz*a))
    
var=1000000
energie=linspace(E0-t0-6*t,E0-t0+6*t,100)

cmpt=linspace(0,100,100)
densite=[]
for i in range(100):
    cmpt[i]=0 
 

for m in range(var):
    k_x=random.uniform(-pi/a,pi/a)
    k_y=random.uniform(-pi/a,pi/a)
    k_z=random.uniform(-pi/a,pi/a)
    E=E0-t0-2*t*(cos(k_x*a)+cos(k_y*a)+cos(k_z*a))
    #for i in range(99):
    #    if energie[i]< E and E < energie[i+1]:
    i=floor((E-(E0-t0-6*t))/(12*t/100.))
    cmpt[i]+=1
    
 

       
#Autre méthode : grille uniforme en 3D
  
"""
"""
for n_x in range(N):
    k_x=-pi/a+n_x*(2*pi/N*a)
    for n_y in range(N):
        k_y=-pi/a+n_y*(2*pi/N*a)
        for n_z in range(N):
            k_z=-pi/a+n_z*(2*pi/N*a)
            E=E0-t0-2*t*(cos(k_x*a)+cos(k_y*a)+cos(k_z*a))
            i=floor((E-(E0-t0-6*t))/(12*t/100.))
            cmpt[i]+=1

                    
#On en déduit directement le nombre d'états d'énergie entre E_i et E_(i+1)
for i in range(100):
    densite.append(cmpt[i]/(12*t/100.))

            
plot(energie,densite)  
"""

    
    

