# -*- coding: utf-8 -*-
"""
Created on Mon Dec 07 17:40:03 2015

@author: Robert
"""



"""
Ce code permet de calculer directement la correction à l'énergie d'un niveau électronique donné (pour un réseau 2D),
grâce à un algorithme de Metropolis. La chaîne de Markov construite est dans l'espace des (l2,p2,r,r').
L'invariance par translation dans la somme a donc été prise en compte.
"""


import numpy 
import scipy
from scipy.integrate import quad, dblquad, tplquad
import math
from math import * 
import random

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


#Estimation de la correction à l'énergie dûe au terme de Hartree-Fock pour des orbitales atomiques gaussiennes

e2=2.3*math.pow(10,-28) 
N=5
E0=13
a=math.pow(10,-10)
d=a/5.

#On choisit n'importe quel couple (l1,p1)
l1=0
p1=0

#Forme de la partie radiale de lorbitale atomique commune à tous les sites
    
def chi_gaussienne_3D(x,y,z):
    return exp(-(x**2+y**2+z**2)/((d/a)**2))
    
#OU une lorentzienne :
    
def chi_lorentzienne_3D(x,y,z):
    return (d/a)**2/(x**2+y**2+z**2+(d/a)**2)
    
x=numpy.linspace(0,5*a,100) 

#chi_1=[(chi_gaussienne(x[i]),chi_lorentzienne(x[i])) for i in range(100)]
#chi_1=[chi_gaussienne(x[i]) for i in range(100)]
#chi_1=[chi_lorentzienne(x[i]) for i in range(100)]


def distance_reseau_2D(x,y,z,x1,y1,z1):
    #§Minimum entre 4 éléments
    var1=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    var2=sqrt((x-x1)**2+(y-(y1-N*a))**2+(z-z1)**2)
    var3=sqrt((x-(x1-N*a))**2+(y-(y1-N*a))**2+(z-z1)**2)
    return min(var1,var2,var3,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 


k_x=numpy.linspace(-pi/a,pi/a,N)
#k_x[O] sera donc égal à -pi/a et k_x[N/2] à  0. On peut translater les indices
k_y=numpy.linspace(-pi/a,pi/a,N)

def terme_facteur_phase(n_x,n_y,l2,p2):
    res=0
    for j_x in range(N):
        for j_y in range(N):
            if (j_x!=n_x) or (j_y!=n_y):
                #rajouter la condition : j est un état occupé de même spin ? 
                #Faire une fonction qui dit si l'état est occupé ou non ?
                res+= -cos((k_x[n_x]-k_x[j_x])*(l1-l2)*a+(k_y[n_y]-k_y[j_y])*(p1-p2)*a)
    #Cas où tous les états sont occupés
    return N**2+res

    
print(terme_facteur_phase(0,0,0,0))

  
#Fonction densité voulue, avec gaussiennes centrées
#Division par N**4 pour normer car somme sur l1,p1,l2,p2
def rho(x0,x1,x2,y0,y1,y2):
    return (chi_gaussienne_3D(x0,x1,x2)*chi_gaussienne_3D(y0,y1,y2))**2/N**2
    
def F(n_x,n_y,l2,p2,x0,x1,x2,y0,y1,y2):
    return terme_facteur_phase(n_x,n_y,l2,p2)/distance_reseau_2D(x0+l1,x1+p1,x2,y0+l2,y1+p2,y2)
    
    
def h(u):
    return min(1,u)

#Taux d'acceptation dans l'algorithme de Metropolis : on prend un taux de transition T(x->y) symétrique
def A(X,Y):
    #return h(rho(Y)*T(Y,X)/(rho(X)*T(X,Y)))
    #return h(rho(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7])*T(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7],X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7])/rho(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7])*T(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6],X[7],Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],Y[7]))
    return h(rho(Y[2],Y[3],Y[4],Y[5],Y[6],Y[7])/rho(X[2],X[3],X[4],X[5],X[6],X[7]))
   
#Taux de transition choisi : à adapter le mieux possible !
def T(n,X,Y):
    return 0

#Taille caractéristique de la boule où l'on tire r_(n+1) à partir de r_n
d1=0.1*d

#Conditions initiales : point de départ dans la chaîne de MArkov :

p1_0=1
p2_0=1
x_0=0.
x_1=0.
x_2=0.
y_0=0.
y_1=0.
y_2=0.

#Nombre de tirages aléatoires dans la marche aléatoire
#Pas vraiment contrôlé ? 
nb=10000.


X=[[0 for m in range(8)] for i in range(10000)]
#Valeurs de la chaîne de Markov


def Delta_Fock_w_s_i_Metropolis(n_x,n_y):
    res=0
    X[0]=[p1_0,p2_0,x_0,x_1,x_2,y_0,y_1,y_2]
    for i in range(9999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires (l,m,r,r')
    
        #On choisit Y selon la loi T(X_i,y)
        Y=[0 for m in range(8)]
        
        var1=random.uniform(0,1)
        if (var1<0.333):
            Y[0]=X[i][0]-1
            if (Y[0]<0):
                Y[0]+=N
                
        if (var1>0.333) and (var1<0.666):
            Y[0]=X[i][0]
        if (var1>0.666):
            Y[0]=X[i][0]+1   
            
        var2=random.uniform(0,1)
        if (var2<0.333):
            Y[1]=X[i][1]-1
            if (Y[1]<0):
                Y[1]+=N
        if (var2>0.333) and (var2<0.666):
            Y[1]=X[i][1]
        if (var2>0.666):
            Y[1]=X[i][1]+1      
            
            
        #Choix d'un nouveau vecteur r (première variable d'intégration) : uiformément choisi dans un petit cube autour du vecteur r précédent
        Y[2]=random.uniform(X[i][2]-(d1/a),X[i][2]+(d1/a))
        Y[3]=random.uniform(X[i][3]-(d1/a),X[i][3]+(d1/a))
        Y[4]=random.uniform(X[i][4]-(d1/a),X[i][4]+(d1/a))
        #Choix d'un nouveau vecteur r' (deuxième variable d'intégration) : uniformément choisi dans un petit cube autour du vecteur r précédent
        Y[5]=random.uniform(X[i][5]-(d1/a),X[i][5]+(d1/a))
        Y[6]=random.uniform(X[i][6]-(d1/a),X[i][6]+(d1/a))
        Y[7]=random.uniform(X[i][7]-(d1/a),X[i][7]+(d1/a))
        
        #Soit U_(i+1) choisi uniformément dans [0,1] (indépendamment du passé)
        nombre=random.uniform(0,1)
        if (nombre<A(X[i],Y)):
            #ceci arrive avec probabilité A(X_i,Y) : taux de transition : dans ce cas on accepte ce Y comme point suivant 
            #de la chaîne de Markov
            X[i+1]=Y
        else:
            X[i+1]=X[i]
        #X[i+1]=[l_(i+1),m_(i+1),r_(i+1),r'_(i+1)]
        #print(X[i+1])
        res+= F(n_x,n_y,X[i+1][0],X[i+1][1],X[i+1][2],X[i+1][3],X[i+1][4],X[i+1][5],X[i+1][6],X[i+1][7])
   
    return (res/nb)
#Cette correction est à mltiplier par e2*a**5
    
#print(Delta_Fock_w_s_i_Metropolis(1))

tab=[[0 for i in range(N)] for j in range(N)]
for m in range(N):
    for q in range(N):
        tab[m][q]=Delta_Fock_w_s_i_Metropolis(m,q)
        print("Delta_Fock_w_s_i_(k_x[{0}]={1},k_y[{2}]={3}) : {4} ".format(m,-pi/a+m*(2*pi/(N*a)),q,-pi/a+q*(2*pi/(N*a)),tab[m][q]))

numpy.savetxt('Correction_energie_2D_bis.txt',tab,newline='\n')


#Graphe en 3D de la correction de l'énergie par le terme de Hartree Fock
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D


#plt.show()

k_x = np.linspace(-pi/a,pi/a, N)
k_x_1=np.outer(k_x,np.ones(np.size(k_x)))

k_y = np.linspace(-pi/a,pi/a,N)
k_y_1=np.outer(np.ones(np.size(k_y)),a*k_y)


correction_energie=[[tab[m][q] for m in range (0,N)] for q in range(0,N)]


ax.plot_surface(k_x_1, k_y_1, correction_energie,rstride=10, cstride=10, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
#La fonction plot_surface trace à partir de 3 matrices A, B et C, l'ensemble des points 
#de coordonnées (A[i][j], B[i][j], C[i][j]) et les relie pour former une surface.
plt.show()
ax.contour(k_x_1, k_y_1, correction_energie,zdir='z')
 

#ax.set_zlim(E0-t0-4*t, E0-t0+4*t)

ax.zaxis.set_major_locator(plt.LinearLocator(10))
ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=10)

plt.show()

"""
X,Y= meshgrid(k_x,k_y)
ax.plot_surface(X,Y,energie_2D)
show()
pcolor(X,Y,energie_2D)# Projection 3D (vue du dessus ; couleurs selon l'intensité)
show()
"""

#save txt numpy : argument array numpy
#lire : load txt
