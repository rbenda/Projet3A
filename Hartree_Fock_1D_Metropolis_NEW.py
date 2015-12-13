# -*- coding: utf-8 -*-
"""
Created on Sat Dec 05 19:17:54 2015

@author: Robert
"""


"""
Ce code permet de calculer directement la correction à l'énergie d'un niveau électronique donné (pour un réseau 1D),
grâce à un algorithme de Metropolis. La chaîne de Markov construite est dans l'espace des (m,r,r').
L'invariance par translation l/m a été prise en compte.
"""

get_ipython().magic(u'pylab inline')
import numpy 
import scipy
from scipy.integrate import quad, dblquad, tplquad
import math
from math import * 
import random
import time

#Estimation de la correction à l'énergie dûe au terme de Hartree-Fock pour des orbitales atomiques gaussiennes

e2=2.3*math.pow(10,-28) 
N=20
E0=13
t0=0.5
t=2
a=math.pow(10,-10)
d=a/20.

#Invariance par translation : on peut fixer la valeur de l à ce que l'on veut
l=0

#Forme de la partie radiale de lorbitale atomique commune à tous les sites
    
def chi_gaussienne_3D(x,y,z):
    return exp(-(x**2+y**2+z**2)/((d/a)**2))
    
#OU une lorentzienne :
    
x0=a/5.
    
def chi_lorentzienne_3D(x,y,z):
    return (d/a)**2/(x**2+y**2+z**2+(d/a)**2)
    
x=numpy.linspace(0,5*a,100) 

#chi_1=[(chi_gaussienne_3D(x[i]),chi_lorentzienne_3D(x[i])) for i in range(100)]
#chi_1=[chi_gaussienne(x[i]) for i in range(100)]
#chi_1=[chi_lorentzienne(x[i]) for i in range(100)]

"""
def integrande_reseau_1D_g(l,m,v_x,v_y,v_z,w_x,w_y,w_z):
    if (w_x-v_x == l-m) and (w_y==v_y) and (w_z==v_z):
        return 0
    else:
        return (chi_gaussienne_3D(w_x,w_y,w_z)*chi_gaussienne_3D(v_x,v_y,v_z))**2*(1/sqrt((w_x-v_x+l-m)**2+(w_y-v_y)**2+(w_z-v_z)**2))
    
def integrande_reseau_1D_l(l,m,v_x,v_y,v_z,w_x,w_y,w_z):
    return (chi_lorentzienne_3D(w_x,w_y,w_z)*chi_lorentzienne_3D(v_x,v_y,v_z))**2*(1/sqrt((w_x-v_x+l-m)**2+(w_y-v_y)**2+(w_z-v_z)**2))
"""
  

def distance_reseau_1D(x,y,z,x1,y1,z1):
    var=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    return min(var,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 


k=numpy.linspace(-pi/a,pi/a,N)
print(k)
#k[O] sera donc égal à -pi/a et k[N/2] à  0. On peut translater les indices

def terme_facteur_phase(n,m):
    res=0
    for j in range(N):
        if (j!=n):
            #rajouter la condition : j est un état occupé de même spin ? Faire une fonction qui dit si l'état est occupé ou non ?
            #print(cos((k[n]-k[j])*(m-l)*a))
            res+= -cos((k[n]-k[j])*(m-l)*a)
    #Cas où tous les états sont occupés (cas général : N_occ - Delta(k_n,occ))
    return N+res

    
print(terme_facteur_phase(0,1))

  
#Fonction densité voulue, avec gaussiennes centrées
#Division par N**2 pour normer car somme sur l et m
def rho(x0,x1,x2,y0,y1,y2):
    return (chi_gaussienne_3D(x0,x1,x2)*chi_gaussienne_3D(y0,y1,y2))**2/N
    
def F(n,m,x0,x1,x2,y0,y1,y2):
    return terme_facteur_phase(n,m)/distance_reseau_1D(x0+l,x1,x2,y0+m,y1,y2)
    
    
def h(u):
    return min(1,u)

#Taux d'acceptation dans l'algorithme de Metropolis : on prend un taux de transition T(x->y) symétrique
def A(X,Y):
    #X et Y sont des tableux à 8 composantes
    #return h(rho(Y)*T(Y,X)/(rho(X)*T(X,Y)))
    #return h(rho(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6])*T(n,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6],X[0],X[1],X[2],X[3],X[4],X[5],X[6])/rho(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6])*T(n,X[0],X[1],X[2],X[3],X[4],X[5],X[6],Y[0],Y[1],Y[2],Y[3],Y[4],Y[5],Y[6]))
    return h(rho(Y[1],Y[2],Y[3],Y[4],Y[5],Y[6])/rho(X[1],X[2],X[3],X[4],X[5],X[6]))
   
#Taux de transition choisi : à adapter le mieux possible !
def T(n,X,Y):
    return 0

#Taille caractéristique de la boule où l'on tire r_(n+1) à partir de r_n
d1=0.001*a
#Plus d1 est grand, plus la zone explorée est grande dans un laps de temps donné



#Conditions initiales : point de départ dans la chaîne de Markov :
#On remarque une grande sensibiité des résultats de correction d'énergie calculée avec cette condition initiale

m0=1
x_0=0.01
x_1=0.
x_2=0.
y_0=0.
y_1=0.
y_2=0.

#Nombre de tirages aléatoires dans la marche aléatoire
#Pas vraiment contrôlé ? 
nb=10000.

#1000000 : insuffisant pour le calcul d'une correction.

X=[[0 for m in range(7)] for i in range(10000)]
#Valeurs de la chaîne de Markov


def Delta_Fock_w_s_i_Metropolis(n):
    res=0
    X[0]=[m0,x_0,x_1,x_2,y_0,y_1,y_2]
    #print(X[0])
    res+= F(n,X[0][0],X[0][1],X[0][2],X[0][3],X[0][4],X[0][5],X[0][6])
    print(res)
    for i in range(9999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires (l,m,r,r')
    
        #On choisit Y selon la loi T(X_i,y)
        Y=[0 for m in range(7)]
        
        #m (deuxième indice de sommation) passe à m+1 avec proba 1/2 et à m-1 avec proba 1/2.
        var2=random.uniform(0,1)
        if (var2<0.333):
            Y[0]=X[i][0]-1
            if (Y[0]<0):
                Y[0]+=N
            
        if (var2>0.333) and (var2<0.666):
            Y[0]=X[i][0]
        if (var2>0.666):
            Y[0]=X[i][0]+1   
            
        #Choix d'un nouveau vecteur r (première variable d'intégration) : uiformément choisi dans un petit cube autour du vecteur r précédent
        Y[1]=random.uniform(X[i][1]-(d1/a),X[i][1]+(d1/a))
        Y[2]=random.uniform(X[i][2]-(d1/a),X[i][2]+(d1/a))
        Y[3]=random.uniform(X[i][3]-(d1/a),X[i][3]+(d1/a))
        #Choix d'un nouveau vecteur r' (deuxième variable d'intégration) : uniformément choisi dans un petit cube autour du vecteur r précédent
        Y[4]=random.uniform(X[i][4]-(d1/a),X[i][4]+(d1/a))
        Y[5]=random.uniform(X[i][5]-(d1/a),X[i][5]+(d1/a))
        Y[6]=random.uniform(X[i][6]-(d1/a),X[i][6]+(d1/a))
        #print(Y)
        #Soit U_(i+1) choisi uniformément dans [0,1] (indépendamment du passé)
        nombre=random.uniform(0,1)
        #print(rho(Y[1],Y[2],Y[3],Y[4],Y[5],Y[6]))
        #print(rho(X[i][1],X[i][2],X[i][3],X[i][4],X[i][5],X[i][6]))
        
        if (nombre<A(X[i],Y)):
            #ceci arrive avec probabilité A(X_i,Y) : taux de transition : dans ce cas on accepte ce Y comme point suivant 
            #de la chaîne de Markov
            X[i+1]=Y
        else:
            X[i+1]=X[i]
        #X[i+1]=[l_(i+1),m_(i+1),r_(i+1),r'_(i+1)]
        #print(X[i+1])
        #print(terme_facteur_phase(n,X[i+1][0]))
        res+= F(n,X[i+1][0],X[i+1][1],X[i+1][2],X[i+1][3],X[i+1][4],X[i+1][5],X[i+1][6])
        #print(res)
    return (res/nb)
#Cette correction est à multiplier par e2*a**5
    
    

n=5
tps1=time.clock()
var =Delta_Fock_w_s_i_Metropolis(n)
tps2=time.clock()
 

print("Delta_Fock_w_s_i_({0}) : {1}  Temps calcul : {2} ".format(n,var,tps2-tps1))   



tab=[0 for i in range(N)]
for m in range(N):
    tps3=time.clock()
    tab[m]=Delta_Fock_w_s_i_Metropolis(m)
    tps4=time.clock()
    print("Delta_Fock_w_s_i_(k[{0}]={1}) : {2} Temps calcul : {3}".format(m,k[m],tab[m],tps4-tps3))


numpy.savetxt('Correction_energie_1D_N=20_nb=1000000_essai1_07_12.txt',tab,newline='\n')
plot(k,tab)
xlabel("k")
ylabel("Delta_HF(k)")


#energie_corrigee=[(E0-t0-2*t*cos(-pi/a+m*(2*pi/(N*a)))+tab[i]) for i in range(N)]
#plot(k,energie_corrigee)


"""
correction_energie=[9.956834077464483457e+00,1.060511466501776212e+01,9.124225664101098587e+00,1.756331837450649402e+01,8.831179766938802800e+00,9.289461306160548659e+00,1.113356831375137723e+01,3.229128590413729682e+00,2.601291496998913377e+00,1.053887144398758480e+01,1.930682997274413282e+01,1.053656427481304192e+01,1.766667668890292830e+01,2.017065349642503591e+01,2.016454664841041122e+01,7.336344856566507922e+00,2.845719310092396093e+00,4.743560815336092418e+00,2.383792207864921409e+00,2.250529804152507651e+01]

energie=[(E0-t0-2*t*cos(-pi/a+m*(2*pi/(N*a)))+correction_energie[i]) for i in range(N)]

plot(k,energie)
"""


#save txt numpy : argument array numpy
#lire : load txt
