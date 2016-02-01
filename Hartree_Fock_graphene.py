# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 21:30:41 2016

@author: Robert
"""


"""
Ce code permet de calculer la correction à l'énergie d'un niveau électronique donné pour le graphène.
"""

get_ipython().magic(u'pylab inline')
import scipy
from scipy.integrate import quad, dblquad
import math
from math import *
from cmath import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy
from numpy import *
import random

N=100
E0=13
t0=0.5
t=2.
a=math.pow(10,-10)
d=a/5.
 

#Construction de k_x et k_y id est de la discrétisation en "k" ; de l'espace réciproque.

#Nombre de tirages aléatoires dans la Zone de Brillouin
NB_tirages=100

for nb in range(NB_tirages):
    #Construction de k_x et k_y  
    alpha=random.uniform(0.,1.)
    beta=0
    #Recherche d'un bon entier "q"
    for q in range(0,floor((4*pi/(3*a))*m)):
        if (q/m-(n/m)*alpha<=1) and (q/m-(n/m)*alpha>=0):
            beta=q/m-(n/m)*alpha
    k_x=(3*a/2)*(alpha+beta)
    k_y=(sqrt(3)*a/2)*(beta-alpha)  
    

def energie_graphene_liant(kx,ky):
    return E0-t0-t*sqrt(3 + 2*(2*cos((3*a/2)*kx)*cos((sqrt(3)*a/2)*ky)+cos(sqrt(3)*a*ky)))
def energie_graphene_antiliant(kx,ky):
    return E0-t0+t*sqrt(3 + 2*(2*cos((3*a/2)*kx)*cos((sqrt(3)*a/2)*ky)+cos(sqrt(3)*a*ky)))
    
var=100000
energie=linspace(E0-t0-3*t,E0-t0+3*t,500)

cmpt=linspace(0,500,500)
densite=[]
for i in range(500):
    cmpt[i]=0 
 
#La zone de Brillouin du graphène est un hexagone inclus dans [-pi/a;pi/a]^2
#Il faut faire des tirages aléatoires de (kx,ky) uniquement dans la 1 Z.B. pour obtenir le bon résultat de densité
#Sinon des points k différents et pas dans la même cellule unité pourraient être tirer et contribuer à la même énergie


#On travaille dans la base (a1*,a2*) du réseau de Bravais réciproque :
A_1=[2*pi/3,-2*pi/sqrt(3)]
A_2=[2*pi/3,2*pi/sqrt(3)]

#alpha=random.uniform(0.,1.)
#beta=random.uniform(0.,1.)
#k1=[(alpha*A_1[i]+beta*A_2[i]) for i in range(2)]


def distance_reseau_2D(x,y,z,x1,y1,z1):
    #Comment périodiser pour le graphène ?
    var1=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    var2=sqrt((x-x1)**2+(y-(y1-N*a))**2+(z-z1)**2)
    var3=sqrt((x-(x1-N*a))**2+(y-(y1-N*a))**2+(z-z1)**2)
    return min(var1,var2,var3,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 



def F(alpha1,alpha2,p1,p2,r1,r2,theta1,theta2,phi1,phi2):
    #return (pi**2)*(2*pi)**2*(pi*(d/a)**2/2)*r1**2*r2**2*sin(theta1)*sin(theta2)/sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*l1*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+l1**2+2*p1*(r1*sin(theta1)*sin(phi1)-r2*sin(theta2)*sin(phi2))+p1**2)
    return r1**2*r2**2*sin(theta1)*sin(theta2)/distance_reseau_2D(r1*sin(theta1)*cos(phi1)+(3*a/2)*(p1-alpha1+p2-alpha2),r1*sin(theta1)*sin(phi1)+(sqrt(3)*a/2)*(p2-alpha2-p1+alpha1),r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2))


#Nombre de tirages aléatoires pour estimer les intégrales
nb=100

def I_Graphene_1(alpha1,alpha2,p1,p2):
    
    Y=[0 for i in range(6)]    
    res=0

    for i in range(99):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        Y[1]=numpy.random.normal(0,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(alpha1,alpha2,p1,p2,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return (pi**2/d)*(a/d)**3*(res/nb)
    
    
def I_Graphene_2(alpha1,alpha2,p1,p2):
    
    Y=[0 for i in range(6)]    
    res=0

    for i in range(99):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(1,(d/a)/2.)
        #Tirage d'une gaussienne centrée en a ? id est en a/a=1 en variables réduites
        Y[1]=numpy.random.normal(1,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(alpha1,alpha2,p1,p2,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return (pi**2/d)*(a/d)**3*(res/nb)
    

def I_Graphene_1_2(alpha1,alpha2,p1,p2):
    
    Y=[0 for i in range(6)]    
    res=0

    for i in range(99):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        #Tirage d'une gaussienne centrée en a ? id est en a/a=1 en variables réduites
        Y[1]=numpy.random.normal(1,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(alpha1,alpha2,p1,p2,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return (pi**2/d)*(a/d)**3*(res/nb)
    
#coefficient de renormalisation
C=1/sqrt(2*N)

def z(kx,ky):
    return 1-2*sin((3*a/2)*kx)*sin((sqrt(3)*a/2)*ky)+j*2*cos((3*a/2)*kx)*sin((sqrt(3)*a/2)*ky)

def z_RL(alpha,beta):
    #Définition de z avec les coordonnées du réseau réciproque
    return 1+cos(2*pi*alpha)+cos(2*pi*beta)+j*(sin(2*pi*alpha)+sin(2*pi*beta))
    
    
def lambda_B_liant(kx,ky):
    return t*z(kx,ky)/(E0-t0-energie_graphene_liant(kx,ky))
    
    
def lambda_B_antiliant(kx,ky):
    return (energie_graphene_antiliant(kx,ky)-(E0-t0))/(t*z(kx,ky).conjugate())
    
  
  
  
def terme_facteur_phase_1(n_x,n_y,alpha1,alpha2,p1,p2):
    res=0
    for j_x in range(N1):
        for j_y in range(N2):
            #Pour les états occupés (k_j_x,k_j_y), avec un électron de même spin que l'état considéré
            phase_x=(k_x[n_x]-k_x[j_x])*(3*a/2)*(p1-alpha1+p2-alpha2)
            phase_y= (k_y[n_y]-k_y[j_y])*(sqrt(3)*a/2)*(alpha1-p1+p2-alpha2)            
            res+= (cos(phase_x)+j*sin(phase_y))
            #Pas de facteurs avec lambda car lambda_A=1
    return res  
    
    
def terme_facteur_phase_2(n_x,n_y,alpha1,alpha2,p1,p2):
    res=0
    for j_x in range(N1):
        for j_y in range(N2):
            #Pour les états occupés (k_j_x,k_j_y), avec un électron de même spin que l'état considéré
            phase_x=(k_x[n_x]-k_x[j_x])*(3*a/2)*(p1-alpha1+p2-alpha2)
            phase_y= (k_y[n_y]-k_y[j_y])*(sqrt(3)*a/2)*(alpha1-p1+p2-alpha2)            
            res+= cos(phase_x)+j*sin(phase_y)
    return res
    
def terme_facteur_phase_1_2(n_x,n_y,alpha1,alpha2,p1,p2):
    res=0
    for j_x in range(N1):
        for j_y in range(N2):
            #Pour les états occupés (k_j_x,k_j_y), avec un électron de même spin que l'état considéré
            phase_x=(k_x[n_x]-k_x[j_x])*(3*a/2)*(p1-alpha1+p2-alpha2)
            phase_y= (k_y[n_y]-k_y[j_y])*(sqrt(3)*a/2)*(alpha1-p1+p2-alpha2)            
            res+= lambda_B_liant(k_x[n_x],k_y[n_y]).conjugate()*lambda_B_liant(k_x[j_x],k_y[j_y])*(cos(phase_x)+j*sin(phase_y))
    return res
    
def terme_facteur_phase_2_1(n_x,n_y,alpha1,alpha2,p1,p2):
    res=0
    for j_x in range(N1):
        for j_y in range(N2):
            #Pour les états occupés (k_j_x,k_j_y), avec un électron de même spin que l'état considéré
            phase_x=(k_x[n_x]-k_x[j_x])*(3*a/2)*(p1-alpha1+p2-alpha2)
            phase_y= (k_y[n_y]-k_y[j_y])*(sqrt(3)*a/2)*(alpha1-p1+p2-alpha2)             
            res +=lambda_B_liant(k_x[n_x],k_y[n_y])*lambda_B_liant(k_x[j_x],k_y[j_y]).conjugate()*(cos(phase_x)+j*sin(phase_y))
    return res
    
    
k_x=linspace(-pi/a,pi/a,10)
k_y=linspace(-pi/a,pi/a,10)
    
#Il reste à faire la somme de ces facteurs de phase sur (alpha1,alpha2) et (p1,p2) e pondérant
#par les intégrales précédentes.
    
#Sommer alpha1, alpha2, p1, p2 d'où à où ? Conditions aux limites ?
    
    
#Vecteur de repliement du feuillet de graphène : C_h=n*a1+m*a2
#n et m définissent le type de nanotube (zigzag, armchair, chiral), et donc aussi son caractère métallique/isolant
C_repliement=[10,10]
    
n=C_repliement[0]
m=C_repliement[1]    

#Tableaux des valeurs des intégrales :
valeurs_I_1=[[[[0 for i in range(m)] for j in range(n)] for k in range(m)] for l in range(n)]
valeurs_I_2=[[[[0 for i in range(m)] for j in range(n)] for k in range(m)] for l in range(n)]
valeurs_I_1_2=[[[[0 for i in range(m)] for j in range(n)] for k in range(m)] for l in range(n)]
valeurs_I_2_1=[[[[0 for i in range(m)] for j in range(n)] for k in range(m)] for l in range(n)]




for p1 in range(n):
    for p2 in range(m):
        for alpha1 in range(n):
            for alpha2 in range(m):
                valeurs_I_1[alpha1][alpha2][p1][p2]=I_Graphene_1(alpha1,alpha2,p1,p2)
                valeurs_I_2[alpha1][alpha2][p1][p2]=I_Graphene_2(alpha1,alpha2,p1,p2)
                valeurs_I_1_2[alpha1][alpha2][p1][p2]=I_Graphene_1_2(alpha1,alpha2,p1,p2)
                valeurs_I_2_1[alpha1][alpha2][p1][p2]=I_Graphene_1_2(alpha1,alpha2,p1,p2)
    
   
   
#null=linspace(0.,0.,N)
#correction_energie=np.outer(null,np.ones(np.size(null)))
    
    
    
for n_x in range():
    for n_y in range():    
        res=0
        for p1 in range(n):
            for p2 in range(m):
                for alpha1 in range(n):
                    for alpha2 in range(m):
                        res+=valeurs_I_1[alpha1][alpha2][p1][p2]*terme_facteur_phase_1(n_x,n_y,alpha1,alpha2,p1,p2)
                        res+=valeurs_I_2[alpha1][alpha2][p1][p2]*terme_facteur_phase_2(n_x,n_y,alpha1,alpha2,p1,p2)

                        res+=valeurs_I_1_2[alpha1][alpha2][p1][p2]*terme_facteur_phase_1_2(n_x,n_y,alpha1,alpha2,p1,p2)
                        res+=valeurs_I_2_1[alpha1][alpha2][p1][p2]*terme_facteur_phase_2_1(n_x,n_y,alpha1,alpha2,p1,p2)
    
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 