# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 12:48:09 2015

@author: Robert
"""
get_ipython().magic(u'pylab inline')
import numpy as np
import scipy
from scipy.integrate import quad, dblquad, tplquad
import math
from math import * 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#Energie pour des électrons libres

h=6.63*math.pow(10,-34)
m_e=9.11*math.pow(10,-31)
e2=2.3*math.pow(10,-28)

def energy_free_electrons(kx,ky,kz): 
    k=sqrt(kx**2+ky**2+kz**2)
    return h**2*k**2/(2*m_e)
    
def energy_free_electrons2(k): 
    return h**2*k**2/(2*m_e)

#Correction à l'énergie trouvée avec le terme de Hartree-Fock

def G(x):
    return 0.5+(1-x**2)/(4*x)*log(abs((1.+x)/(1.-x)))


def energy_Hartree_Fock_electrons(kx,ky,kz,k_F): 
    k=sqrt(kx**2+ky**2+kz**2)
    
    #Methode avec l'intégrale triple
    #def integrande1(r1,theta1,phi1):
     
      #  return (r1**2*sin(theta1))/(k**2+r1**2-2*r1*(sin(theta1)*(kx*cos(phi1)+ky*sin(phi1))+kz*cos(theta1)))
    
    #return energy_free_electrons(kx,ky,kz)-(4*pi*e**2)/(2*pi)**3*tplquad(lambda r1, theta1, phi1 : integrande1(r1,theta1,phi1), 0.,2*pi,lambda theta1 : 0., lambda theta1 : pi, lambda r1 : 0., lambda r1 : k_F)
    
    #Méthode en admettant la formule dans Ashcroft & Mermin
    return energy_free_electrons(kx,ky,kz)-((2*e2)/pi)*k_F*G(k/k_F)
    
def energy_Hartree_Fock_electrons2(k,k_F):   
    return energy_free_electrons2(k)-((2*e2)/pi)*k_F*G(k/k_F)
    
    

k_F=math.pow(10,8)
print(((2*e2)/pi)*k_F)
k_1=np.linspace(0,2*k_F,100)
 
G_1=[G(k_1[i]/k_F) for i in range(100)]
#plot(k_1,G_1)

#Que se passe-t-il au voisinage de k_F ?
k_2=np.linspace(0.7*k_F,1.3*k_F,100)

E_free=[energy_free_electrons2(k_1[i]) for i in range(100)]
E_corrige=[energy_Hartree_Fock_electrons2(k_1[i],k_F) for i in range(100)]
 

#plot(k_1, E_free, k_1, E_corrige)

#hold on?

#plot(k_1,E_free,'g')
#show()
#plot(k_1,E_corrige,'b')


#ax.set_ylim(-3*k_F,-2*k_F)

E_tot=[(E_free[i],E_corrige[i]) for i in range(100)]

plot(k_1,E_tot)
xlabel("k")
ylabel("E(k) according to Hartree-Fock")
show() 
hold()


k0=0.5*math.pow(10,8)

def energie_a_k_fixe(kF):
    return h**2*k0**2/(2*m_e)-((2*e2)/pi)*kF*G(k0/kF)
    
k_F_1=np.linspace(k0,k_F,100)
E_a_k_fixe=[energie_a_k_fixe(k_F_1[i]) for i in range(100)]

plot(k_F_1,E_a_k_fixe,'r')
xlabel("k_F")
ylabel("E(k_F) according to Hartree-Fock pour k=k0=0.5*10^8 fixe ")
show()
hold()


#Derivee de l'énergie de Hartree Fock par rapport à k

k_1=linspace(0.99999999*k_F,1.00000001*k_F,1000)
derivee_energy_Hartree_Fock_electrons=[((energy_Hartree_Fock_electrons2(k_1[i+1],k_F)-energy_Hartree_Fock_electrons2(k_1[i],k_F))/((k_1[i+1]-k_1[i]))) for i in range(0,999)]

#Divergence logarithmique en kF : zoomer 10 fois ne fait que doubler la taille du pic /etc...

plot(k_1[0:999],derivee_energy_Hartree_Fock_electrons)
xlabel("k")
ylabel("dE_corrige(k)/dk selon Hartree-Fock pour kF=10^8 ")
hold()

#Evaluation de Delta(k) par une méthode de Monte-Carlo


def distance(x,y,z,x1,y1,z1):
    return sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)    
    
#Tirage d'une variable uniforme dans la boule de centre k et de rayon kF

#↨Nombre de tirages
var=100000.

def delta_k_Monte_Carlo(kx,ky,kz,kF):
    res=0
    for i in range(10000):
        #Tirages de 3 variables aléatoires indépendantes, uniformes sur [kx-kF,kx+kF], [ky-kF,ky+kF] et [kz-kF,kz+kF]
        X1=random.uniform(kx-kF,kx+kF)
        X2=random.uniform(ky-kF,ky+kF)
        X3=random.uniform(kz-kF,kz+kF)
        if (distance(X1,X2,X3,kx,ky,kz)<=kF) :
            #(X1,X2,X3)= v est dans la boule de centre k et de rayon kF
            res += 1/(X1**2+X2**2+X3**2)
    return ((4*pi*e2)/(2*pi)**3)*(2*kF)**3*(res/var)
    
def energy_Hartree_Fock_Monte_Carlo(kx,ky,kz,kF): 
    return energy_free_electrons(kx,ky,kz)-delta_k_Monte_Carlo(kx,ky,kz,kF)
        

k_F=math.pow(10,8)

k_1=np.linspace(0,2*k_F,100)

E_free=[energy_free_electrons2(k_1[i]) for i in range(100)]
E_corrige=[energy_Hartree_Fock_Monte_Carlo(k_1[i],0,0,k_F) for i in range(100)]
 
E_tot=[(E_free[i],E_corrige[i]) for i in range(100)]

plot(k_1,E_tot)
xlabel("k")
ylabel("E(k) according to Hartree-Fock with Monte-Carlo")
hold()


#Evaluation de Delta(k) en calculant l'intégrale triple : sphériques ou cartésiennes ??
#def delta_k_integrale_triple(kx,ky,kz,kF):
#    k=sqrt(kx**2+ky**2+kz**2)
#    def integrande1(r1,theta1,phi1):
#        return (r1**2*sin(theta1))/(k**2+r1**2-2*r1*(sin(theta1)*(kx*cos(phi1)+ky*sin(phi1))+kz*cos(theta1)))
#    return energy_free_electrons(kx,ky,kz)-(4*pi*e**2)/(2*pi)**3*tplquad(lambda r1, theta1, phi1 : integrande1(r1,theta1,phi1), 0.,2*pi,lambda theta1 : 0., lambda theta1 : pi, lambda r1 : 0., lambda r1 : kF)
  
"""
def delta_k_integrale_triple_bis(kx,ky,kz,kF):
    k=sqrt(kx**2+ky**2+kz**2)
    def integrande_1(r1,theta1,phi1):
        return (r1**2*sin(theta1))/(k**2+r1**2-2*r1*(sin(theta1)*(kx*cos(phi1)+ky*sin(phi1))+kz*cos(theta1)))
        
    def integrande_2(theta1,phi1):
        return quad(lambda r1 : integrande_1(r1,theta1,phi1), 0.,kF)
        
    def integrande_3(phi1):
        return quad(lambda theta1 : integrande_2(theta1,phi1),0.,pi)
        
    return (4*pi*e**2)/(2*pi)**3*quad(lambda phi1 : integrande_3(phi1),0.,2*pi)
   
def energy_Hartree_Fock_integrale_triple_bis(kx,ky,kz,kF): 
    return energy_free_electrons(kx,ky,kz)-delta_k_integrale_triple_bis(kx,ky,kz,kF)
           
k_1=np.linspace(0,2*k_F,100)

E_free=[energy_free_electrons2(k_1[i]) for i in range(100)]
E_corrige=[energy_Hartree_Fock_integrale_triple_bis(k_1[i],0,0,k_F) for i in range(100)]
 
E_tot=[(E_free[i],E_corrige[i]) for i in range(100)]

plot(k_1,E_tot)
"""
 
