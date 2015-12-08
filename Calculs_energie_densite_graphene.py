# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:39:44 2015

@author: Robert
"""


get_ipython().magic(u'pylab inline')


import scipy
from scipy.integrate import quad, dblquad
import math
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D

k_x = np.linspace(-pi/a,pi/a, 200)
k_x_1=np.outer(a*k_x,ones(size(k_x)))

k_y = np.linspace(-pi/a,pi/a,200)
k_y_1=np.outer(ones(size(k_y)),a*k_y)

var_1=np.outer(cos((3*a/2)*k_x), cos((sqrt(3)*a/2)*k_y))
#print(var_1)
var_2=np.outer(ones(size(k_y)), cos(sqrt(3)*a*k_y))
#print(var_2)

def energy_graphene_antiliant(u,v):
    return E0-t0+t*sqrt(3+2*(2*cos((3*a/2)*u)*cos((sqrt(3)*a/2)*v)+cos(sqrt(3)*a*v)))

def energy_graphene_liant(u,v):
    return E0-t0-t*sqrt(3+2*(2*cos((3*a/2)*u)*cos((sqrt(3)*a/2)*v)+cos(sqrt(3)*a*v)))

#Energie pour un feuillet de graphène : énergie des états liants et antiliants
energie_graphene_antiliant=(E0-t0+t*sqrt(3+2*(2*var_1+var_2)))

energie_graphene_antiliant_bis=[[energy_graphene_antiliant(-pi/a+i*(2*pi/a)/200,-pi/a+j*(2*pi/a)/200) for i in range (0,200)] for j in range(0,200)]

energie_graphene_liant=(E0-t0-t*sqrt(3+2*(2*var_1+var_2)))

energie_graphene_liant_bis=[[energy_graphene_liant(-pi/a+i*(2*pi/a)/200,-pi/a+j*(2*pi/a)/200) for i in range (0,200)] for j in range(0,200)]


ax.plot_surface(k_x_1, k_y_1, energie_graphene_antiliant_bis,rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)

#ax.plot_surface(k_x_1, k_y_1, energie_graphene_liant,rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.plot_surface(k_x_1, k_y_1, energie_graphene_liant_bis,rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_zlim(E0-t0-3*t, E0-t0+3*t)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()

X1,Y1= meshgrid(k_x,k_y)
ax.plot_surface(X1,Y1,energie_graphene_liant_bis)
#ax.plot_surface(X1,Y1,energie_graphene_antiliant_bis)
show()

#Visualisation des points de Dirac vus du dessus
pcolor(X1,Y1,energie_graphene_liant)# Projection
show()
