# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 13:36:49 2015

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

"""
import sympy
from sympy import Integral,limit,integrate
from sympy.abc import x
from sympy import Symbol,cos,erf,exp,sqrt


x=Symbol('x')
print(limit(1/x,x,0))

x=Symbol('x')
print(integrate(exp(-x**2)*erf(x), (x, 0, 1)))

#sympy.integrate pour les intégrales impropres
x=Symbol('x')
print(integrate(sqrt(x), (x, 0, 1)))
"""



N=100
E0=13
t0=0.5
t=2
a=math.pow(10,-10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D


#plt.show()

k_x = np.linspace(-pi/a,pi/a, 100)
k_x_1=np.outer(a*k_x,np.ones(np.size(k_x)))

k_y = np.linspace(-pi/a,pi/a,100)
k_y_1=np.outer(np.ones(np.size(k_y)),a*k_y)

x_1=np.outer(np.cos(a*k_x),np.ones(np.size(k_x)))
#print(x_1)
x_2=np.outer(np.ones(np.size(k_y)),np.cos(a*k_y))
#print(x_2)

#Energie pour un réseau carré 2D de pas a
energie_2D=(E0-t0-2*t*(x_1+x_2))
#print(energie_2D)

def energy_2D(u,v):
    return E0-t0-2*t*(cos(a*u)+cos(a*v))

energie_2D_bis=[[energy_2D(-pi/a+i*(2*pi/a)/100,-pi/a+j*(2*pi/a)/100) for i in range (0,100)] for j in range(0,100)]


ax.plot_surface(k_x_1, k_y_1, energie_2D,rstride=10, cstride=10, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
#La fonction plot_surface trace à partir de 3 matrices A, B et C, l'ensemble des points 
#de coordonnées (A[i][j], B[i][j], C[i][j]) et les relie pour former une surface.
#xlabel("k_x (10^(10) m^(-1)")
#ylabel("k_y (10^(10) m^(-1)")
plt.show()
ax.contour(k_x_1, k_y_1, energie_2D,zdir='z')


ax.set_zlim(E0-t0-4*t, E0-t0+4*t)

ax.zaxis.set_major_locator(plt.LinearLocator(10))
ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()

"""
X,Y= meshgrid(k_x,k_y)
#ax.plot_surface(X,Y,energie_2D_bis)
#show()
pcolor(X,Y,energie_2D)# Projection 3D (vue du dessus ; couleurs selon l'intensité)
#xlabel("k_x (m^(-1))")
#ylabel("k_y (m^(-1))")
show()
"""


def integrande(z,E) :
    return 1/sqrt((1-z**2)*(1-((E0-t0-2*t*z-E)/(2*t))**2))

#Fonction test
#def g(x) : 
#    return x*log(x)

#print(quad(lambda x: g(x), 0,1)[0])

# A quoi correspond quad(lambda x: g(x), 0,1)[1] ??


#Densité d'états d'énergie en 2D
def D_2D(E):
    if E > E0-t0 and E<E0-t0+4*t:
        return N*N/(2*t*pi*pi)*quad(lambda z: integrande(z,E), -1, (E0-t0+2*t-E)/(2*t))[0] #Intégrale de -1 à alpha(E)+1
    #(alpha(E)<0)
    if E < E0-t0 and E>E0-t0-4*t:
        return N*N/(2*t*pi*pi)*quad(lambda z: integrande(z,E), (E0-t0-2*t-E)/(2*t),1)[0] #Intégrale de alpha(E)-1 à 1
    #(alpha(E)>0)
    if E < E0-t0-4*t or E>E0-t0+4*t:
        return 0

 
print("D_2D(E0-t0) : {0}".format(D_2D(E0-t0)))

nb=1000.0
E_2D = np.linspace(E0-t0-5*t,E0-t0+5*t,nb)
densites_2D= [D_2D(E0-t0-5*t+i*(10*t)/1000.) for i in range(0,1000,1)]


#Densité d'états d'énergie en 2D
plot(E_2D,densites_2D,label="D(E)_2D")
xlabel("Energy (eV)")
ylabel("DOS D(E)")
show()
hold()

E_2D_plus = np.linspace(E0-t0,E0-t0+5*t,nb)
essai=[N**2*(1/sqrt(E0-t0+i*(5*t)/100.)) for i in range(1000)]
#plot(E_2D_plus,essai)

#print((N/pi)*quad(lambda z: D_2D(E0-t0-t+2*t*cos(z)),0,pi)[0])

#Vérification de la normalisation : l'intégrale totale de D_2D(E) doit valoir N**2.
#Pas de découpe de l'intégrale de Riemann à signaler ? 
print("Intégrale de D_2D jusqu'à {1} : {0}".format(quad(lambda E : D_2D(E),E0-t0-4*t,E0-t0-0.01*t)[0],E0-t0-0.01*t))
print("Intégrale de D_2D jusqu'à {1} : {0}".format(quad(lambda E : D_2D(E),E0-t0+0.01*t,E0-t0+4*t)[0],E0-t0+0.01*t))

#Reprogrammation de l'intégrale de Riemann de l'intégrale entre E0-t0-4*t et E0-t0
pas=1000.
subdivision=[(E0-t0-4*t+i*(4*t/pas)) for i in range(1000)]

integrale_riemann=0 
for i in range(1000-1):
    integrale_riemann += (subdivision[i+1]-subdivision[i])*D_2D(subdivision[i+1])
    
print("Integrale par somme de Riemann : {0}".format(integrale_riemann))
#OK : on trouve bien quelque chose tendant vers 5000=N^2/2 quand le pas de la subdivision tend vers 0 ("pas" tend vers l'infini)

print(quad(lambda z : 1/sqrt(z),0,1)[0])


#x=Symbol('x')
#print("Intégrale de D_2D jusqu'à {1} : {0}".format(integrate(D_2D(x),(x,E0-t0-4*t,E0-t0)),E0-t0))

#x=Symbol('x')
#print(sympy.limit(D_2D(x),x,13.))


def integrande_2(z,E) :
    return np.arccos(((E0-t0-E)/(2*t))-z)/sqrt((1-z**2))

#Nombre d'états en 2D
def N_2D(E):
    if E > E0-t0 and E<=E0-t0+4*t:
        return (N**2/pi**2)*(pi*np.arccos((E0-t0+2*t-E)/(2*t))+quad(lambda z: integrande_2(z,E), -1, (E0-t0+2*t-E)/(2*t))[0])
    if E < E0-t0 and E>=E0-t0-4*t:
         return N**2/(2*pi**2)*quad(lambda z: integrande_2(z,E), (E0-t0-2*t-E)/(2*t),1)[0]
    if (E<E0-t0-4*t):
        return 0
    if (E>E0-t0+4*t):
        return N_2D(E0-t0+4*t)
    
    
nb=1000.0

#Visualisation de la discontinuité de la pente en E0-t0-4t
E_2D = linspace(E0-t0-4.5*t,E0-t0-3.5*t,nb)
nb_etats_2D= [N_2D(E0-t0-4.5*t+i*(1*t)/1000.) for i in range(0,1000,1)]
plot(E_2D,nb_etats_2D)
xlabel("Energy (eV)")
ylabel("N_2D(E)")
show()

E_2D = np.linspace(E0-t0-5*t,E0-t0+5*t,nb)
nb_etats_2D= [N_2D(E0-t0-5*t+i*(10*t)/1000.) for i in range(0,1000,1)]

print("N_2D(E0-t0+4*t) : {0}".format(N_2D(E0-t0+4*t)))
print("N_2D(E0-t0) : {0}".format(N_2D(E0-t0)))
print("N_2D(E0-t0-4*t) : {0}".format(N_2D(E0-t0-4*t)))

plt.plot(E_2D,nb_etats_2D)
xlabel("Energy (eV)")
ylabel("N_2D(E)")
#Même résultat en traçant la dérivée de D_2D(E)?

essai=[10*N*(sqrt(E0-t0+i*(5*t)/100.)) for i in range(1000)]
#plot(E_2D_plus,essai)


#print(N_2D(E0-t0+4*t-(8*t)/1000000.))

#print(N_2D(E0-t0-4*t+(1)*(8*t/nb))-N_2D(E0-t0-4*t+0*(8*t/nb)))



nb=100.
densite_2D_bis=[(N_2D(E0-t0-4*t+(i+2)*(8*t/nb))-N_2D(E0-t0-4*t+(i+1)*(8*t/nb))) for i in range(0,98,1)]
E_2D_bis=linspace(E0-t0-4*t+(8*t/nb),E0-t0+4*t,98)

plot(E_2D_bis,densite_2D_bis)



