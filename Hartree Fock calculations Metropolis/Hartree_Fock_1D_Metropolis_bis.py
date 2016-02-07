# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:30:00 2015

@author: Robert
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 22:59:27 2015

@author: Robert
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 10:19:55 2015

@author: Robert
"""


"""
Ce code permet de calculer chaque intégrale I_1D(m) (et chaque facteur de phase Delta(m)), intervenant dans la correction de l'énergie pour un réseau 1D,
grâce à un algorithme de Metropolis : la chaîne de Markov construite à (m) fixé est dans l'espace des (r,r') à 6 dimensions

Permet de tracer le graphe des I_1D(m)*Delta(n,m) à n fixé, pour tous les m=1..N. (N chaînes de MArkov construites)
On retrouve ainsi la décroissance en 1/m à partir du site 0.
"""

get_ipython().magic(u'pylab inline')
import numpy 
import scipy
from scipy.integrate import quad, dblquad, tplquad
import math
from math import * 
import random

#Estimation de la correction à l'énergie dûe au terme de Hartree-Fock pour des orbitales atomiques gaussiennes

e2=2.3*math.pow(10,-28) 
N=50
E0=13
a=math.pow(10,-10)
d=a/4.
l=0

#Forme de la partie radiale de lorbitale atomique commune à tous les sites
    
def chi_gaussienne_3D(x,y,z):
    return exp(-(x**2+y**2+z**2)/((d/a)**2))
    
#OU une lorentzienne :
    
x0=a/5.
    
def chi_lorentzienne_3D(x,y,z):
    return (d/a)**2/(x**2+y**2+z**2+(d/a)**2)
    
x=numpy.linspace(0,5*a,100) 

#chi_1=[(chi_gaussienne(x[i]),chi_lorentzienne(x[i])) for i in range(100)]
#chi_1=[chi_gaussienne(x[i]) for i in range(100)]
#chi_1=[chi_lorentzienne(x[i]) for i in range(100)]

"""
def integrande_reseau_1D_g(m,v_x,v_y,v_z,w_x,w_y,w_z):
    if (w_x-v_x == l-m) and (w_y==v_y) and (w_z==v_z):
        return 1
    else:
        return (chi_gaussienne_3D(w_x,w_y,w_z)*chi_gaussienne_3D(v_x,v_y,v_z))**2*(1/sqrt((w_x-v_x+l-m)**2+(w_y-v_y)**2+(w_z-v_z)**2))
    
def integrande_reseau_1D_l(l,m,v_x,v_y,v_z,w_x,w_y,w_z):
    return (chi_lorentzienne_3D(w_x,w_y,w_z)*chi_lorentzienne_3D(v_x,v_y,v_z))**2*(1/sqrt((w_x-v_x+l-m)**2+(w_y-v_y)**2+(w_z-v_z)**2))
"""
    
  
def distance_reseau_1D(x,y,z,x1,y1,z1):
    var=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    return min(var,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 




NB=100
#Remplissage en électrons

k=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]


def terme_facteur_phase(n,m):
    res=0
    for j in range(NB/2):
        if (j!=n):
            #rajouter la condition : j est un état occupé de même spin ? Faire une fonction qui dit si l'état est occupé ou non ?
            #print(cos((k[n]-k[j])*(m-l)*a))
            res+= -cos((k[n]-k[(N-NB/2)/2+j])*(m-l)*a)
    #Cas où tous les états sont occupés
    return NB+res

    




#Conditions initiales : point de départ dans la chaîne de Markov :
x_0=0.
x_1=0.
x_2=0.
y_0=0.
y_1=0.
y_2=0.


#Fonction densité voulue : sera la mesure de invariante de la chaîne de MArkov construite
def rho(m,x0,x1,x2,y0,y1,y2):
    #return integrande_reseau_1D_g(l,m,x0,x1,x2,y0,y1,y2)
    return (chi_gaussienne_3D(x0,x1,x2)*chi_gaussienne_3D(y0,y1,y2))**2
    
    
#Multiplier juste par le facteurde phase pour le MEtropolis plus général ?    
def F(m,x0,x1,x2,y0,y1,y2):
    return 1/distance_reseau_1D(x0+l,x1,x2,y0+m,y1,y2)
    
def h(u):
    return min(1,u)

#Taux d'acceptation dans l'algorithme de Metropolis
def A(m,X,Y):
    #X et Y sont des tableux à 8 composantes
    #return h(rho(Y)*T(Y,X)/(rho(X)*T(X,Y)))
    #return h(rho(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7])*T(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7],X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7])/rho(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7])*T(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7]))
    return h(rho(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])/rho(m,X[0],X[1],X[2],X[3],X[4],X[5]))
   
#Taux de transition choisi : à adapter le mieux possible !
def T(n,X,Y):
    return 0
 

  
#Nombre de tirages aléatoires dans la marche aléatoire
#Pas vraiment contrôlé ? 
nb=10000.

d1=0.1*a

X=[[0 for m in range(6)] for i in range(10000)]
#Valeurs de la chaîne de Markov

def I_1D_Metropolis(m):
    res=0
    X[0]=[x_0,x_1,x_2,y_0,y_1,y_2]
    for i in range(9999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires (r,r') 
        #On choisit Y selon la loi T(X_i,y)
        Y=[0 for m in range(6)]
        #Choix d'un nouveau vecteur r (première variable d'intégration) : uniformément choisi dans un petit cube autour du vecteur r précédent
        Y[0]=random.uniform(X[i][0]-(d1/a),X[i][0]+(d1/a))
        Y[1]=random.uniform(X[i][1]-(d1/a),X[i][1]+(d1/a))
        Y[2]=random.uniform(X[i][2]-(d1/a),X[i][2]+(d1/a))
        
        #Choix d'un nouveau vecteur r' (deuxième variable d'intégration) : uniformément choisi dans un petit cube autour du vecteur r précédent
        Y[3]=random.uniform(X[i][3]-(d1/a),X[i][3]+(d1/a))
        Y[4]=random.uniform(X[i][4]-(d1/a),X[i][4]+(d1/a))
        Y[5]=random.uniform(X[i][5]-(d1/a),X[i][5]+(d1/a))
        
        #Soit U_(i+1) choisi uniformément dans [0,1] (indépendamment du passé)
        nombre=random.uniform(0,1)
        #print(rho(l,m,X[i][0],X[i][1],X[i][2],X[i][3],X[i][4],X[i][5]))
        if (nombre<A(m,X[i],Y)):
            #ceci arrive avec probabilité A(X_i,Y) : taux de transition : dans ce cas on accepte ce Y comme point suivant 
            #de la chaîne de Markov
            X[i+1]=Y
        else:
            X[i+1]=X[i]
        #X[i+1]=[r_(i+1),r'_(i+1)]
        #print(X[i+1])
        if (i>=5000):
            #On ne fait la moyenne des termes le long de la chaîne de Markov qu'une fois 
            #la "thermalisation" effectuée, i.e. que la mesure stationnaire est presque atteinte.
            res+= F(m,X[i+1][0],X[i+1][1],X[i+1][2],X[i+1][3],X[i+1][4],X[i+1][5])
   
    return (res/4999.)
    
    
#Ce qui compte :X[nb]=(r_nb,r'_nb) : où en est arrivé la chaîne de Markov après suffisamment d'étapes.
 
  
n=0
tab=[0 for i in range(N)]
for m_1 in range(N):
    tab[m_1]=I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1)
    print("I_1D({0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)))
    print("Delta({2},{0}) : {1} ".format(m_1,terme_facteur_phase(n,m_1),n))
    print("I_1D({0})*Delta({2},{0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1),n))

M=[i for i in range(N)]
plot(M,tab)


"""
def Delta_Fock_w_s_i(n):   
    res=0
    for l_1 in range (0,N,1):
        for m_1 in range (0,N,1):
            #print("terme_facteur_phase{1} {2} :{0}".format(terme_facteur_phase(n,l_1,m_1),l_1,m_1))
            #print("I_1D_MonteCarlo:{0}".format(I_1D_MonteCarlo[l_1][m_1]))
            res+= I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1)
            #-(e2/N**2)*terme_facteur_phase(n,l_1,m_1)
    return res



for m in range(N):
    print("Delta_Fock_w_s_i_(k[{0}]={1}) : {2} ".format(m,k[m],Delta_Fock_w_s_i(m)))
"""

#save txt numpy : argument array numpy
#lire : load txt