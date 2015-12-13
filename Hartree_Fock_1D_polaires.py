# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 18:25:55 2015

@author: Robert
"""

get_ipython().magic(u'pylab inline')
import numpy 
import scipy
import math
from math import * 
import random
import time


e2=2.3*math.pow(10,-28) 
N=50
E0=13
a=math.pow(10,-10)

#d : inférieur ou égal à a/4. pour une gaussienne (ainsi 1 % de recouvrement des gaussiennes
# de 2 orbitales atomiques localisée voisines en a/2)
d=a/4.



#Cette fonction sert en coordonnées cartésiennes uniquement
#Distance sur un "anneau" 1D périodisé (site 0 = site N)
def distance_reseau_1D(x,y,z,x1,y1,z1):
    var=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    return min(var,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 


#Nombre d'ELECTRONS mis dans le système
#On remplit les états de façon non magnétique : en partant de l'état fondamental ; 
#2 électrons (spins up et down) par état
NB=N
#NB=N_occ/2 où N-occ est le nombre d'ETATS occupés

#Demi-remplissage : NB=N 
#(N niveaux en tout ; 2 électrons par niveaux : donc 2N électrons au total au maximum)

k=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]
#k=numpy.linspace(-(2*pi/(N*a))*floor(N/2.),(2*pi/(N*a))*floor(N/2.),N)
#print(k)
#k doit contenir la valeur 0 : k(floor(N/2.))=0 nécessairement

def terme_facteur_phase(n,m):
    res=0
    for j in range(NB/2):
    #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
        if (j!=n):
            #rajouter la condition : j est un état occupé de même spin ? Faire une fonction qui dit si l'état est occupé ou non ?
            #print(cos((k[n]-k[j])*(m-l)*a))
            res+= -cos((k[n]-k[N/4+j])*m*a)
    #Cas où tous les états sont occupés
    return NB+res

    
def F(m,r1,r2,theta1,theta2,phi1,phi2):
    return (pi**2)*(2*pi)**2*(pi*(d/a)**2/2)*r1**2*r2**2*sin(theta1)*sin(theta2)/sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*m*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+m**2)
    
#Facteur exp(-2*r1**2/(d/a)**2)*exp(-2*r1**2/(d/a)**2) dans l'intégrande : contenu dans
#le tirage de r1 et r2 selon des gaussiennes    
#Premier facteur : facteur de renormalisation de la densité de la gaussienne
    
nb=100000.

r1_0=0.
r2_0=0.
theta1_0=0.
theta2_0=0.
phi1_0=0.
phi2_0=0.


def I_1D(m):
    
    Y=[0 for m in range(6)]    
    res=0
    Y=[r1_0,r2_0,theta1_0,theta2_0,phi1_0,phi2_0]
    for i in range(99999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        Y[1]=numpy.random.normal(0,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return (res/nb)
#Facteur K**4 * a**5 (où K est la partie angulaire de l'orbitale atomique)


"""
n=N/2
tab=[0 for i in range(N)]
for m_1 in range(N):
    tab[m_1]=I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1)
    #print("I_1D({0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)))
    #print("Delta({2},{0}) : {1} ".format(m_1,terme_facteur_phase(n,m_1),n))
    print("I_1D({0})*Delta({2},{0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1),n))

M=[i for i in range(N)]
plot(M,tab)
"""



#Avec cette méthode, on recalcule pour chaque m 
#I_1D(m) avec nb tirages du 6-uplet (rho1,rho2,theta1,theta2,phi1,phi2)
def Delta_Hartree_Fock(n):
    res2=0
    tab_I_1D=[0 for i in range(N)]
    for m in range(N):
        tab_I_1D[m]=I_1D(m)
        res2+= tab_I_1D[m]*terme_facteur_phase(n,m)
        #Le calcul de I_1D(m) nécessite à chaque m nb tirages...
        
    abscisses_m=[i for i in range(N)]
    plot(abscisses_m,tab_I_1D)
    return (e2/N)*res2
    
print(Delta_Hartree_Fock(2))


nb_tirages=10000.

#Avec cette méthode, on évalue pour chaque tirage la valeur de l'intégrande de I_1D
#pour chaque m.
def Delta_Hartree_Fock_bis(n):
    res=[0 for i in range(N)]
    
    Y=[0 for i in range(6)]

    for i in range(9999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        Y[1]=numpy.random.normal(0,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(-pi,pi) 
        Y[5]=random.uniform(-pi,pi)
     
        tab=[0 for i in range(N)]
        
        for m in range(N):
            res[m]+=F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
            tab[m] = F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
            
        #Vérification que pour un tirage des variables donné, m -> F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]) évolue bien en 1/|m|
        #abscisses_m=[i for i in range(10)]
        #print(Y)
        #print(F(1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]))
        #plot(abscisses_m,tab[0:10])
        
    I_1D_estime=[0 for i in range(N)]
    I_1D_estime[m]=res[m]/nb_tirages
    
    resultat_correction_energie=0
    for k in range(N):
        resultat_correction_energie += I_1D_estime[k]*terme_facteur_phase(n,k)
        
    return (e2/N)*resultat_correction_energie


#print(Delta_Hartree_Fock_bis(0))


"""
correction_energie=[0 for i in range(N)]
   
for n in range(N):
    tps1=time.clock()
    correction_energie[n]=Delta_Hartree_Fock_bis(n)
    tps2=time.clock()
    print("Delta_Hartree_Fock({0})= {1}  Temps calcul : {2}".format(n,correction_energie[n],tps2-tps1))

plot(k,correction_energie)
xlabel("k")
ylabel("Delta_HF(k)")

numpy.savetxt('Correction_energie_1D_polaires_N=50_nb=50000_NB=50_essai1_13_12.txt',tab,newline='\n')
"""