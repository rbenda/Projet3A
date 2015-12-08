# -*- coding: utf-8 -*-
"""
Éditeur de Spyder
 
Ceci est un script temporaire.
"""

# coding: utf-8

# In[1]:

get_ipython().magic(u'pylab inline')
import math
def g(x):
    return math.pow(x,2)

print(g(2))
print(g(10.01))

a1=[(0+i*0.1) for i in range(0,10,1)]
a1.append(2)
print(a1)

print(sum([0,1]))

print(3**3)##a**n = a*a*...*a n fois


# In[3]:

from math import pi
import numpy
N=100
E0=13
t0=0.5
t=2
a=math.pow(10,-10)

# Energie pour le réseau d'atomes en 1D avec un pas de a
energie_1D=[(E0-t0-2*t*cos((2*pi*i/(N*a))*a)) for i in range(-N/2,N/2,1)]
abscisses=[(2*pi*i/(N*a)) for i in range(-N/2,N/2,1)]

plot(abscisses,energie_1D)
legend(['E(k)'])


# In[8]:

E = linspace(E0-t0-2*t,E0-t0+2*t,100)#Dernier argument : nombre de points découpant l'intervalle.

#Densité d'états en 1D
z = N/(2*t)*(1/sqrt((1-((E-(E0-t0))/(2*t))**2)))

plt.plot(E, z,label= "D(E)")
xlabel("Energie (eV)")

#print(numpy.arccos((E-(E0-t0))/(2*t)))

nb_etats_1D=N-(N/pi)*numpy.arccos((E-(E0-t0))/(2*t)) # Nombre d'états d'énergie inférieure ou égale à E

plot(E,nb_etats_1D,label="N(E)")


# In[9]:

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

x = 10 * np.outer(np.cos(u), np.sin(v)) #Matrice : produit tensoriel des deux vecteurs cos(u) et sin(v)
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
#ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='b')

#plt.show()

k_x = np.linspace(-pi/a,pi/a, 100)
k_x_1=np.outer(a*k_x,ones(size(k_x)))

k_y = np.linspace(-pi/a,pi/a,100)
k_y_1=np.outer(ones(size(k_y)),a*k_y)

x_1=np.outer(cos(a*k_x),ones(size(k_x)))
#print(x_1)
x_2=np.outer(ones(size(k_y)),cos(a*k_y))
#print(x_2)

#Energie pour un réseau carré 2D de pas a
energie_2D=(E0-t0-2*t*(x_1+x_2))
#print(energie_2D)

def energy_2D(u,v):
    return E0-t0-2*t*(cos(a*u)+cos(a*v))

energie_2D_bis=[[energy_2D(-pi/a+i*(2*pi/a)/100,-pi/a+j*(2*pi/a)/100) for i in range (0,100)] for j in range(0,100)]


ax.plot_surface(k_x_1, k_y_1, energie_2D_bis,rstride=10, cstride=10, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#La fonction plot_surface trace à partir de 3 matrices A, B et C, l'ensemble des points 
#de coordonnées (A[i][j], B[i][j], C[i][j]) et les relie pour former une surface.
plt.show()
ax.contour(k_x_1, k_y_1, energie_2D,zdir='z')


ax.set_zlim(E0-t0-4*t, E0-t0+4*t)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()


X,Y= meshgrid(k_x,k_y)
ax.plot_surface(X,Y,energie_2D)
show()
pcolor(X,Y,energie_2D)# Projection 3D (vue du dessus ; couleurs selon l'intensité)
show()


# In[27]:


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


# In[11]:

import scipy
from scipy.integrate import quad

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
    else:
        return 0

nb=1000.0
E_2D = linspace(E0-t0-5*t,E0-t0+5*t,nb)
densites_2D= [D_2D(E0-t0-5*t+i*(10*t)/1000.) for i in range(0,1000,1)]

#Densité d'états d'énergie en 2D
plot(E_2D,densites_2D,label="D(E)_2D")


print((N/pi)*quad(lambda z: D_2D(E0-t0-t+2*t*cos(z)),0,pi)[0])


# In[30]:


def integrande_2(z,E) :
    return numpy.arccos(((E0-t0-E)/(2*t))-z)/sqrt((1-z**2))

#Nombre d'états en 2D
def N_2D(E):
    if E > E0-t0 and E<=E0-t0+4*t:
        return (N**2/pi**2)*(pi*numpy.arccos((E0-t0+2*t-E)/(2*t))+quad(lambda z: integrande_2(z,E), -1, (E0-t0+2*t-E)/(2*t))[0])
    if E < E0-t0 and E>E0-t0-4*t:
         return N**2/(2*pi**2)*quad(lambda z: integrande_2(z,E), (E0-t0-2*t-E)/(2*t),1)[0]
    if (E<E0-t0-4*t):
        return 0
    if (E>E0-t0+4*t):
        return N_2D(E0-t0+4*t)
    
    
nb=1000.0

#Visualisation de la discontinuité de la pente en E0-t0-4t
#E_2D = linspace(E0-t0-4.5*t,E0-t0-3.5*t,nb)
#nb_etats_2D= [N_2D(E0-t0-4.5*t+i*(1*t)/1000.) for i in range(0,1000,1)]

E_2D = linspace(E0-t0-5*t,E0-t0+5*t,nb)
nb_etats_2D= [N_2D(E0-t0-5*t+i*(10*t)/1000.) for i in range(0,1000,1)]

plot(E_2D,nb_etats_2D)
#Même résultat en traçant la dérivée de D_2D(E)?

#print(N_2D(E0-t0+4*t-(8*t)/1000000.))

#print(N_2D(E0-t0-4*t+(1)*(8*t/nb))-N_2D(E0-t0-4*t+0*(8*t/nb)))

"""
nb=100.
densite_2D_bis=[(N_2D(E0-t0-4*t+(i+2)*(8*t/nb))-N_2D(E0-t0-4*t+(i+1)*(8*t/nb))) for i in range(0,98,1)]
E_2D_bis=linspace(E0-t0-4*t+(8*t/nb),E0-t0+4*t,98)

plot(E_2D_bis,densite_2D_bis)
"""


# In[32]:

def alpha(E):
    return (E0-t0-E)/(2*t)

#Nombre d'états en 3D
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


E_3D = linspace(E0-t0-7*t,E0-t0+7*t,50)
nb_etats_3D= [N_3D(E0-t0-7*t+i*(14*t)/nb2) for i in range(0,50,1)]

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
        #return (N/pi)*(quad(lambda z: D_2D(E+2*t*cos(z)),0,pi))[0]

"""densite_etats_3D=[(N_3D(E0-t0-6*t+(i+1)*(12*t)/nb2)-N_3D(E0-t0-6*t+i*(12*t)/nb2))/(12*t/nb2) for i in range(0,nb2-1,1)]
"""

E_3D = linspace(E0-t0-6*t,E0-t0+6*t,50)
densite_etats_3D= [D_3D(E0-t0-6*t+i*(12*t)/50.) for i in range(0,50,1) ]

plot(E_3D,densite_etats_3D)


# In[55]:

from scipy.integrate import quad
def integrand(x, a, b):
    return 1/sqrt(x)
a = 2
b = 3
I = quad(integrand, 0, 1,args=(a,b))
print(I)


# In[ ]:

#Estimation de la correction à l'énergie dûe au terme de Hartree-Fock pour des orbitales atomiques gaussiennes


