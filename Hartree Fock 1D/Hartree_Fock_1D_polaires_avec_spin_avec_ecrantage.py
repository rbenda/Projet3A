# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 15:53:20 2016

@author: Robert
"""


"""
Ce code permet de calculer la correction à l'énergie d'un niveau électronique donné (pour un réseau 1D),
en écrivant les intégrales I_1D(m) en coordonnées sphériques et en tirant les variables radiales
selon des lois gaussiennes et les variables angulaires uniformément.
L'invariance par translation dans la somme a été prise en compte.
Erreur d'indice (j au lieu de (N-NB/2)/2+j) corrigée.
Prise en compte du remplissage en spin (et donc de possibles remplissages magnétiques) dans
le terme de facteur de phase. Ceci permet également de distinguer les cas N pairs des cas N impairs,
ce qui n'était pas le cas dans les programmes précédents.

Ce code prend en compte l'écrantage en remplacant le potentiel de Coulomb par un potentiel de Yukawa.
"""


get_ipython().magic(u'pylab inline')
import numpy 
import scipy
import math
from math import * 
import random
import time


e2=2.3*math.pow(10,-28) 
N=200
E0=13.
t0=0.5
t=2.
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
NB=102
#NB=2*N_occ où N_occ est le nombre d'ETATS occupés


#Remplissage des N états d'énergie avec les NB électrons, en tenant compte du spin.
#La somme du nombre de cases valant 1 dans tab_spin_up et dans tab_spin_down doit valoir NB.
tab_spin_up=[0 for i in range(N)]
tab_spin_down=[0 for i in range(N)]

if (NB%2==0):
    if ((NB/2)%2==1):
    #Cas où NB=4l+2
        for k in range((N-NB/2)/2,(N-NB/2)/2+NB/2):
            tab_spin_up[k]=1
            tab_spin_down[k]=1
             
    if ((NB/2)%2==0):
    #Cas où NB=4l
        #Remplissage des états centraux de plus basse énergie et totalement remplis (avec 4l-2 électrons)
        for k in range((N-NB/2)/2+1,(N-NB/2)/2+NB/2):
            tab_spin_up[k]=1
            #if (k>=(N-NB/2)/2+NB/8) and (k<=(N-NB/2)/2+3*NB/8):
            tab_spin_down[k]=1
        #Rempissage ferromagnétique: (deux électrons restants)
        tab_spin_up[(N-NB/2)/2]=1
        tab_spin_up[(N-NB/2)/2+NB/2]=1
        
        """
        #Rempissage anti-ferromagnétique:
        tab_spin_down[(N-NB/2)/2]=1
        tab_spin_down[(N-NB/2)/2+NB/2]=1

        #Autre remplissage:
        tab_spin_up[(N-NB/2)/2]=1
        tab_spin_down[(N-NB/2)/2+NB/2]=1
        """
if (NB%2==1):
    if (((NB-1)/2)%2==0):
    #Cas où NB=4l+1
        #Remplissage des états centraux de plus basse énergie et totalement remplis (avec 4l-2 électrons)
        for k in range((N-NB/2)/2+1,(N-NB/2)/2+(NB-3)/2+1):
            tab_spin_up[k]=1
            tab_spin_down[k]=1
        #Choix d'un remplissage pour les 3 électrons résiduels sur les deux états (symétriques)
        #de plus haute énergie
        tab_spin_down[(N-NB/2)/2]=1
        tab_spin_up[(N-NB/2)/2]=1
        tab_spin_up[(N-NB/2)/2+(NB-3)/2+1]=1
        
    if (((NB-1)/2)%2==1):
    #Cas où NB=4l+3
        #Remplissage des états centraux de plus basse énergie et totalement remplis (avec 4l+2 électrons)
        for k in range((N-NB/2)/2,(N-NB/2)/2+(NB-1)/2):
            tab_spin_up[k]=1
            tab_spin_down[k]=1 
        #Choix de l'état de l'électron qui reste (deux niveaux d'énergie symétriques possibles, et deux spins possibles)
        tab_spin_up[(N-NB/2)/2-1]=1
        #OU
        #tab_spin_down[(N-NB/2)/2-1]=1
        #OU
        #tab_spin_up[(N-NB/2)/2+(NB-1)/2]=1
        #OU        
        #tab_spin_down[(N-NB/2)/2+(NB-1)/2]=1
        
        

#Demi-remplissage : NB=N 
#(N niveaux en tout ; 2 électrons par niveaux : donc 2N électrons au total au maximum)

k=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]
#k=numpy.linspace(-(2*pi/(N*a))*floor(N/2.),(2*pi/(N*a))*floor(N/2.),N)
#print(k)
#k doit contenir la valeur 0 : k(floor(N/2.))=0 nécessairement

def terme_facteur_phase(n,sigma_n,m):
    
#tab_spin_up,tab_spin_down : à passer en paramètres de cette fonction ??
    
    res=0 
    
    if (sigma_n==1):
    #Par convention : 1 désigne up
        res=0
        #Somme des termes de Fock sans interaction et d'Hartree sans interaction
        for j in range(N):
            #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
            if (j!=n):
                res+= -tab_spin_up[j]*cos((k[n]-k[j])*m*a)
                #tab_spin_up[j] vaudra 0 pour des états non occupés, 
                #ou pour des états occupés par des spins "down" seuls ; et 1 sinon
       
        return NB-tab_spin_up[n]+res
        
        
    if (sigma_n==-1):
    #Par convention : -1 désigne spin down

        res=0
        #Somme des termes de Fock sans interaction et d'Hartree sans interaction
        for j in range(N):
            #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
            if (j!=n):
                res+= -tab_spin_down[j]*cos((k[n]-k[j])*m*a)
                #tab_spin_down[j] vaudra 0 pour des états non occupés, 
                #ou pour des états occupés par des spins "up" seuls ; et 1 sinon
      
        return NB-tab_spin_down[n]+res
            
print(terme_facteur_phase(1,-1,1))

N1=200.
NB1=102.
kF=(pi/(2*a))*(NB1/N1)
#Energie de Fermi en tight-binding :
E_F=E0-t0-2*t*cos(kF*a)


lambda_ecrantage=[0 for i in range(5)]

#Valeur initiale : densitée calculée en tight-bonding 1D
lambda_1=(pi/sqrt(e2))*sqrt(2*t*a*abs(sin(kF*a)))/kF
#Pas homogène 
#Attention : e2 = q2/4*pi*epsilon0 !!
lambda_courant=lambda_1

for etape in range(5):
    #Calcul auto-cohérent tenant compte de l'écrantage pour effacer la singularité au niveau de Fermi

    #Longueur caractéristique d'écrantage de Thomas-Fermi, déduite de la densité d'états au niveau de Fermi
    #de façon auto-cohérente (actualisée à chaque étape)

    
    lambda_ecrantage[etape]=lambda_courant
    print("Lambda_ecrantage Etape {0} : {1}".format(etape,lambda_courant))
        
    def F(m,r1,r2,theta1,theta2,phi1,phi2):
        #var=sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*m*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+m**2)
        #distance=min(var,sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*(m-N)*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+(m-N)**2))
      
        return r1**2*r2**2*sin(theta1)*sin(theta2)*exp(-(1/lambda_ecrantage[etape])*distance_reseau_1D(r1*sin(theta1)*cos(phi1)-m,r1*sin(theta1)*sin(phi1),r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2)))/distance_reseau_1D(r1*sin(theta1)*cos(phi1)-m,r1*sin(theta1)*sin(phi1),r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2))
            
            
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
   
    k_F=k[(N+NB/2)/2]
    print("k_F={0}".format(k_F))
   
    correction_energie_spin_up=[0 for i in range(N)]
    correction_energie_spin_down=[0 for i in range(N)]
                    
    valeurs_I_1D=[0 for i in range(N)]
    #On calcule les valeurs de I_1D une fois pour toutes (en faisant beaucoup de tirages pour un m donné et ainsi avoir peu d'incertitudes)
    #(elles sont indépendantes de n ; n n'intervient que dans le facteur de phase)
    for m in range(N):
        valeurs_I_1D[m]=I_1D(m)
      
    for n in range(N):
        tps1=time.clock()
                            
        res_spin_up=0
        res_spin_down=0
                            
        for m in range(N):
            #print("terme_facteur_phase({0},{1}) = {2}".format(v,u,terme_facteur_phase(v,u)))
            res_spin_up += valeurs_I_1D[m]*terme_facteur_phase(n,1,m)
            res_spin_down += valeurs_I_1D[m]*terme_facteur_phase(n,-1,m)
            
            correction_energie_spin_up[n]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res_spin_up
            correction_energie_spin_down[n]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res_spin_down
            
            tps2=time.clock()
            #print("Delta_Hartree_Fock({0})= {1}  Temps calcul : {2}".format(v,correction_energie[v],tps2-tps1))
            
            
    correction_energie=[(correction_energie_spin_up[i],correction_energie_spin_down[i]) for i in range(N)]
            
    plot(k,correction_energie)
    xlabel("k")
    ylabel("Delta_HF(k), polarisee en spin en eV")
    show()
    hold()
                                
                            
    #Energie calculee en tight-binding : E0, t0 et t sont exprimés en eV
    energie_sans_correction=[(E0-t0-2*t*cos(k[i]*a)) for i in range(N) ]
                                
    energie_corrigee_spin_up=[(E0-t0-2*t*cos(k[i]*a)+correction_energie_spin_up[i]) for i in range(N)]
                                
    comparaison_energie_up=[(energie_sans_correction[i],energie_corrigee_spin_up[i]) for i in range(N)]
    
    energie_corrigee_spin_down=[(E0-t0-2*t*cos(k[i]*a)+correction_energie_spin_down[i]) for i in range(N)]
                                
    comparaison_energie_down=[(energie_sans_correction[i],energie_corrigee_spin_down[i]) for i in range(N)]
                                
                                
    """
    niveau_Fermi_non_corrige=[(E0-t0-2*t*cos(k_F*a)) for i in range(N)]
                                
    niveau_Fermi_corrige=[energie_corrigee[(N-NB/2)/2] for i in range(N)]
                                
    largeur_bande_etats_occupes_apres_correction=energie_corrigee[(N+NB/2)/2]-energie_corrigee[N/2]
    Delta_occ=largeur_bande_etats_occupes_apres_correction
                                
    largeur_bande_etats_vides_apres_correction=energie_corrigee[0]-energie_corrigee[(N-NB/2)/2]
    Delta_vide=largeur_bande_etats_vides_apres_correction
                                
                                
    Delta_occ_1=energie_sans_correction[(N+NB/2)/2]-energie_sans_correction[N/2]
    Delta_vide_1=energie_sans_correction[0]-energie_sans_correction[(N-NB/2)/2]
    """
                                
    plot(k,comparaison_energie_up) 
    #text(-10000000000, 17, r'$\Delta L_{occ}, \Delta L_{empty}$',fontsize=17)
    xlabel("k")
    ylabel("E(k), E(k)_corrigee pour spins UP en eV")
    show()
    hold()
    #scatter([k_F],[0])
                                
                                
                                
    plot(k,comparaison_energie_down) 
    #text(-10000000000, 17, r'$\Delta L_{occ}, \Delta L_{empty}$',fontsize=17)
    xlabel("k")
    ylabel("E(k), E(k)_corrigee pour spins DOWN en eV")
    show()
    hold()
                                
    comparaison_energie_up_down=[(energie_corrigee_spin_up[i],energie_corrigee_spin_down[i]) for i in range(N)]
                                
                                
    plot(k,comparaison_energie_up_down) 
    #text(-10000000000, 17, r'$\Delta L_{occ}, \Delta L_{empty}$',fontsize=17)
    xlabel("k")
    ylabel("E(k)_corrigee pour UP, E(k)_corrigee pour DOWN, en eV")
    show()
    hold()
    """
    #Largeur de bandes des états occupés, largeur de bande des états vides :
    print("Largeurs de bande des états vides et occupés après correction par Hartree-Fock :")
    print("Post-correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ,Delta_vide))
    print("Post-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ/Delta_vide))
    print("Avant correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ_1,Delta_vide_1))
    print("Avant-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ_1/Delta_vide_1))
    """
    #print("Ces rapport doivent évoluer de façon monotone avec le remplissage NB, à N fixé")
                                    
                                    
    #Derivée de l'énergie corrigée
    #energie_corrigee=[(E0-t0-2*t*cos(k[i]*a)+correction_energie[i]) for i in range(1000)]
    derivee_energie_corrigee=[((energie_corrigee_spin_up[i+1]-energie_corrigee_spin_up[i])/(k[i+1]-k[i])) for i in range(0,N-1)]
    #plot(k[(N+NB/2)/2-10:(N+NB/2)/2+10],derivee_energie_corrigee)
    plot(k[0:N-1],derivee_energie_corrigee)
    #title("Screening included")
    xlabel("k")
    ylabel("dE_corr.(k)/dk en eV.m")
    show()
    hold()
                                
                                
    #Actualisation de la longueur caractéristique d'écrantage de Thomas-Fermi grâce au calcul de
    #la densité d'états au niveau de Fermi pour la structure de bande corrigée calculée par Hartree-Fock
    derivee_energie_kF=derivee_energie_corrigee[(N+NB/2)/2]
                                
    lambda_courant=(pi/sqrt(e2))*sqrt(abs(derivee_energie_kF))/kF
                                
                                