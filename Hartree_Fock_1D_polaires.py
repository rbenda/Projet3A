# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 18:25:55 2015

@author: Robert
"""

"""
Ce code permet de calculer la correction à l'énergie d'un niveau électronique donné (pour un réseau 1D),
en écrivant les intégrales I_1D(m) en coordonnées sphériques et en tirant les variables radiales
selon des lois gaussiennes et les variables angulaires uniformément.
L'invariance par translation dans la somme a été prise en compte.
Erreur d'indice (j au lieu de (N-NB/2)/2+j) corrigée
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
t0=0.5
t=2
a=math.pow(10,-10)

#d : inférieur ou égal à a/4. pour une gaussienne (ainsi 1 % de recouvrement des gaussiennes
# de 2 orbitales atomiques localisée voisines en a/2)
d=a/4.

#TESTER d'AUTRES VALEURS DE d et voir comment varie l'ordre de grandeur de la correction de Fock


#Distance sur un "anneau" 1D périodisé (site 0 = site N)
def distance_reseau_1D(x,y,z,x1,y1,z1):
    var=sqrt((x-(x1-N))**2+(y-y1)**2+(z-z1)**2)
    return min(var,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 


#Nombre d'ELECTRONS mis dans le système
#On remplit les états de façon non magnétique : en partant de l'état fondamental ; 
#2 électrons (spins up et down) par état
NB=8
#NB=2*N_occ où N_occ est le nombre d'ETATS occupés

#Demi-remplissage : NB=N 
#(N niveaux en tout ; 2 électrons par niveaux : donc 2N électrons au total au maximum)

k=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]
#k=numpy.linspace(-(2*pi/(N*a))*floor(N/2.),(2*pi/(N*a))*floor(N/2.),N)
#print(k)
#k doit contenir la valeur 0 : k(floor(N/2.))=0 nécessairement

def terme_facteur_phase(n,m):
    res=0
    """
    #Somme des termes d'auto-interaction de Fock et d'Hartree uniquement
    for j in range(NB/2):
    #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
        if ((N-NB/2)/2+j==n):
            #Terme d'auto-interaction de Fock
            res+= -cos((k[n]-k[(N-NB/2)/2+j])*m*a)
    if (n < (N-NB/2)/2) or (n >= (N+NB/2)/2):
        #Dans ce cas l'état "k_n" n'est pas occupé
        return res
    if (n >= (N-NB/2)/2) or (n < (N+NB/2)/2):
        #Dans ce cas l'état "k_n" est occupé
        #Terme d'Hartree sans-autointeraction : (N-1) ; terme d'auto-interaction de Hartree 1, si k_n est un état occupé
        return 1+res  
    #return res
    """
    
    #Somme des termes de Fock sans interaction et d'Hartree sans interaction
    for j in range(NB/2):
    #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
        if ((N-NB/2)/2+j!=n):
            res+= -cos((k[n]-k[(N-NB/2)/2+j])*m*a)
    if (n < (N-NB/2)/2) or (n >= (N+NB/2)/2):
        #Dans ce cas l'état "k_n" n'est pas occupé
        return NB+res
    if (n >= (N-NB/2)/2) or (n < (N+NB/2)/2):
        #Dans ce cas l'état "k_n" est occupé
        return NB-1+res
    return res
    
    
def F(m,r1,r2,theta1,theta2,phi1,phi2):
    #var=sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*m*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+m**2)
    #distance=min(var,sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*(m-N)*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+(m-N)**2))
    #return r1**2*r2**2*sin(theta1)*sin(theta2)/distance
    return r1**2*r2**2*sin(theta1)*sin(theta2)/distance_reseau_1D(r1*sin(theta1)*cos(phi1)-m,r1*sin(theta1)*sin(phi1),r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2))
    
    
#Facteur exp(-2*r1**2/(d/a)**2)*exp(-2*r1**2/(d/a)**2) dans l'intégrande : contenu dans
#le tirage de r1 et r2 selon des gaussiennes    
#Premier facteur : facteur de renormalisation de la densité de la gaussienne
    
nb=1000.


def I_1D(m):
    
    Y=[0 for i in range(6)]    
    res=0

    for i in range(999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(0,(d/a)/sqrt(2.))
        Y[1]=numpy.random.normal(0,(d/a)/sqrt(2.))
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return  (pi**2/d)*(a/d)**3*(res/nb)
#Facteur K**4 * a**5 (où K est la partie angulaire de l'orbitale atomique)


"""

tab=[0 for i in range(N)]
for m_1 in range(N):
    tab[m_1]=I_1D(m_1)
    #print("I_1D({0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)))
    #print("Delta({2},{0}) : {1} ".format(m_1,terme_facteur_phase(n,m_1),n))
    #print("I_1D({0})*Delta({2},{0}) : {1} ".format(m_1,I_1D_Metropolis(m_1)*terme_facteur_phase(n,m_1),n))
    print("I_1D({0}) : {1} ".format(m_1,I_1D(m_1)))


M=[i for i in range(N)]
plot(M,tab)
"""

"""
#Estimation de l'écart type sur I_1D(m) pour un nombre de tirages donné
m=2
  
#Nombre de calculs de la correction (=nb de fois que l'on a simulé une chaîne de Markov)
p=10.

#Estimateur de la moyenne
mu=0.

#Somme des réalisations
sum=0.

#Estimteur de la moyenne du carré
moy=0.

#Somme des carrés des réalisations
sum2=0.

#Estimateur de l'écart type
sigma=0.

for i in range(10):

    tps1=time.clock()
    var =I_1D(m)
    tps2=time.clock()
    sum +=var
    sum2+=var**2
    print("I_1D({0}) : {1}  Temps calcul : {2} ".format(m,var,tps2-tps1))   

mu=sum/p
moy=sum2/p
sigma=sqrt(moy-mu**2)
print("Ecart-type : {0}  Moyenne : {1}".format(sigma,mu))
print("Pourcentage : {0} %".format(100*sigma/mu))
"""


"""
#Avec cette méthode, on recalcule pour chaque m 
#I_1D(m) avec nb tirages du 6-uplet (rho1,rho2,theta1,theta2,phi1,phi2)
def Delta_Hartree_Fock(n):
    res2=0
    tab_I_1D=[0 for i in range(N)]
    for m in range(N):
        tab_I_1D[m]=I_1D(m)*terme_facteur_phase(n,m)
        res2+= tab_I_1D[m]
        #Le calcul de I_1D(m) nécessite à chaque m nb tirages...
        
    #abscisses_m=[i for i in range(N)]
    #plot(abscisses_m,tab_I_1D)
    return (e2/N)*res2
    
#print(Delta_Hartree_Fock(2))


nb_tirages=100000.

#Avec cette méthode, on évalue pour chaque tirage la valeur de l'intégrande de I_1D
#pour chaque m.
#Rajouter variable de spin en argument ?
def Delta_Hartree_Fock_bis(n):
    res=[0 for i in range(N)]
    
    Y=[0 for i in range(6)]

    for i in range(99999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        Y[1]=numpy.random.normal(0,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(-pi,pi) 
        Y[5]=random.uniform(-pi,pi)
     
        #tab=[0 for i in range(N)]
        
        for m in range(N):
            res[m]+=F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
            #tab[m] = F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
            
        #Vérification que pour un tirage des variables donné, m -> F(m,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]) évolue bien en 1/|m|
        #abscisses_m=[i for i in range(10)]
        #print(Y)
        #print(F(1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]))
        #plot(abscisses_m,tab[0:10])
        
    I_1D_estime=[0 for i in range(N)]
    
    for k in range(N):
        I_1D_estime[k]=(4*pi**3/d)*(a/d)**3*(res[k]/nb_tirages)
    
    resultat_correction_energie=0
    
    for k in range(N):
        resultat_correction_energie += I_1D_estime[k]*terme_facteur_phase(n,k)
        
    return (e2/N)*resultat_correction_energie
"""


""" 
n=0
 
#Estimation de l'écart type pour un nombre de tirages donné
    
#Nombre de calculs de la correction (=nb de fois que l'on a simulé une chaîne de Markov)
p=10.

#Estimateur de la moyenne
mu=0.

#Somme des réalisations
sum=0.

#Estimteur de la moyenne du carré
m=0.

#Somme des carrés des réalisations
sum2=0.

#Estimateur de l'écart type
sigma=0.

for i in range(10):

    tps1=time.clock()
    var =Delta_Hartree_Fock_bis(n)
    tps2=time.clock()
    sum +=var
    sum2+=var**2
    print("Delta_Hartree_Fock_bis({0}) : {1}  Temps calcul : {2} ".format(n,var,tps2-tps1))   

mu=sum/p
m=sum2/p
sigma=sqrt(m-mu**2)
print("Ecart-type : {0}  Moyenne : {1}".format(sigma,mu))
print("Pourcentage : {0} %".format(100*sigma/mu))
"""



k_F=k[(N+NB/2)/2]
print("k_F={0}".format(k_F))

correction_energie=[0 for i in range(N)]

valeurs_I_1D=[0 for i in range(N)]
#On calcule les valeurs de I_1D une fois pour toutes (en faisant beaucoup de tirages pour un m donné et ainsi avoir peu d'incertitudes)
#(elles sont indépendantes de n ; n n'intervient que dans le facteur de phase)
for m in range(N):
    valeurs_I_1D[m]=I_1D(m)

"""    
M=[i for i in range(N)]
plot(M,valeurs_I_1D)
show()
hold()
"""
   
#print(valeurs_I_1D)
    
for v in range(N):
    tps1=time.clock()
    
    res=0
    for u in range(N):
        #print("terme_facteur_phase({0},{1}) = {2}".format(v,u,terme_facteur_phase(v,u)))
        res+= valeurs_I_1D[u]*terme_facteur_phase(v,u)
        
    correction_energie[v]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res
    
    tps2=time.clock()
    #print("Delta_Hartree_Fock({0})= {1}  Temps calcul : {2}".format(v,correction_energie[v],tps2-tps1))



plot(k,correction_energie)
xlabel("k")
ylabel("Delta_HF(k) en eV")
show()
hold()

 
#numpy.savetxt('Correction_energie_1D_polaires_N=50_nb=10000_NB=50_20_12_methode_1.txt',correction_energie,newline='\n')


#Energie calculee en tight-binding : E0, t0 et t sont exprimés en eV
energie_sans_correction=[(E0-t0-2*t*cos(k[i]*a)) for i in range(N) ]

#energie_corrigee_ajustee=[(E0-t0-2*t*cos(k[i])+2*t*(correction_energie[i]-correction_energie[N/2])/((correction_energie[0]+correction_energie[1])/2-correction_energie[N/2])) for i in range(N)]

energie_corrigee=[(E0-t0-2*t*cos(k[i]*a)+correction_energie[i]) for i in range(N)]

niveau_Fermi_non_corrige=[(E0-t0-2*t*cos(k_F*a)) for i in range(N)]

niveau_Fermi_corrige=[energie_corrigee[(N-NB/2)/2] for i in range(N)]

energie=[(energie_sans_correction[i],energie_corrigee[i]) for i in range(N)]



largeur_bande_etats_occupes_apres_correction=energie_corrigee[(N+NB/2)/2]-energie_corrigee[N/2]
Delta_occ=largeur_bande_etats_occupes_apres_correction

largeur_bande_etats_vides_apres_correction=energie_corrigee[0]-energie_corrigee[(N-NB/2)/2]
Delta_vide=largeur_bande_etats_vides_apres_correction


Delta_occ_1=energie_sans_correction[(N+NB/2)/2]-energie_sans_correction[N/2]
Delta_vide_1=energie_sans_correction[0]-energie_sans_correction[(N-NB/2)/2]


plot(k,energie) 
#text(-10000000000, 17, r'$\Delta L_{occ}, \Delta L_{empty}$',fontsize=17)
xlabel("k")
ylabel("E(k), E(k)_corrigee en eV")
show()
hold()
#scatter([k_F],[0])

#LArgeur de bandes des états occupés, largeur de bande des états vides :
print("Largeurs de bande des états vides et occupés après correction par Hartree-Fock :")
print("Post-correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ,Delta_vide))
print("Post-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ/Delta_vide))
print("Avant correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ_1,Delta_vide_1))
print("Avant-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ_1/Delta_vide_1))

#print("Ces rapport doivent évoluer de façon monotone avec le remplissage NB, à N fixé")


#Derivée de l'énergie corrigée
#energie_corrigee=[(E0-t0-2*t*cos(k[i]*a)+correction_energie[i]) for i in range(1000)]
derivee_energie_corrigee=[((energie_corrigee[i+1]-energie_corrigee[i])/(k[i+1]-k[i])) for i in range((N+NB/2)/2-5,(N+NB/2)/2+5)]
plot(k[(N+NB/2)/2-5:(N+NB/2)/2+5],derivee_energie_corrigee)
xlabel("k")
ylabel("dE_corr.(k)/dk en eV.m")
#Divergence logarithmique en kF : zoomer 10 fois ne fait que doubler la taille du pic /etc...
