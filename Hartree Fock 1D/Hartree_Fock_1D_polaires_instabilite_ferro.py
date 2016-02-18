# -*- coding: utf-8 -*-
"""
Created on Mon Feb 08 13:18:22 2016

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

Ce code inclut le calcul de l'énergie de différentes configurations magnétiques : E(M) où
M est l'aimantation totale : M=Nb(spins UP)-Nb(spins DOWN) (en remplissant par niveaux croissants d'énergie
à la fois pour les spins UP et pour les spins DOWN).
L'énergie E(M) est la somme des énergies de chaque état occupé (x 2 si les deux types de spin y sont présents),
l'énergie étant bien sûr l'énergie CORRIGEE calculée par Hartree-Fock.
On veut ainsi savoir si la méthode de Hartree-Fock entraîne une instabilité ferromagnétique,
id est une valeur de l'aimantation M non nulle correspondant à un minimum de l'énergie E(M).

"""


get_ipython().magic(u'pylab inline')
import numpy 
import scipy
import math
from math import * 
import random
import time


e2=2.3*math.pow(10,-28) 
N=30
E0=13
t0=0.5
t=2
a=math.pow(10,-10)

#d : inférieur ou égal à a/4. pour une gaussienne (ainsi 1 % de recouvrement des gaussiennes
#de 2 orbitales atomiques localisée voisines en a/2)
d=a/4.

#TESTER d'AUTRES VALEURS DE d et voir comment varie l'ordre de grandeur de la correction de Fock


#Distance sur un "anneau" 1D périodisé (site 0 = site N)
def distance_reseau_1D(x,y,z,x1,y1,z1):
    var=sqrt((x-(x1-N))**2+(y-y1)**2+(z-z1)**2)
    return min(var,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 
    


k_vect=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]
#k doit contenir la valeur 0 : k(floor(N/2.))=0 nécessairement

#Energie calculee en tight-binding : E0, t0 et t sont exprimés en eV
energie_sans_correction=[(E0-t0-2*t*cos(k_vect[i]*a)) for i in range(N)]

def F(m,r1,r2,theta1,theta2,phi1,phi2):
    return r1**2*r2**2*sin(theta1)*sin(theta2)/distance_reseau_1D(r1*sin(theta1)*cos(phi1)-m,r1*sin(theta1)*sin(phi1),r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2))
    
    
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



correction_energie_spin_up=[0 for i in range(N)]
correction_energie_spin_down=[0 for i in range(N)]
        
valeurs_I_1D=[0 for i in range(N)]
#On calcule les valeurs de I_1D une fois pour toutes (en faisant beaucoup de tirages pour un m donné et ainsi avoir peu d'incertitudes)
#(elles sont indépendantes de n ; n n'intervient que dans le facteur de phase)
for m in range(N):
    valeurs_I_1D[m]=I_1D(m)

#NB : nombre d'électrons au total : fixé pour le calcul d'un profil d'énergie M -> E(M)
for NB in range(10,15):
    
    print("Trace de l'energie en fonction de l'aimantation pour un nombre de particules fixé égal à NB={0}".format(NB))
    
    energie_M=[0 for i in range(0,NB+1,1)]
    #On aura energie_M[k] qui correspond à l'énergie totale pour une aimantation -NB+2*k
    
    for NB_DOWN in range(0,NB+1):
        
        NB_UP=NB-NB_DOWN
        #Nombre d'ELECTRONS de spins UP et DOWN mis dans le système : NB_UP et NB_DOWN
        
        #Aimantation M (dans [-N,N])
        M=NB_UP-NB_DOWN
        
        indice_aimantation=(M-(-NB))/2
        #M-(-NB) est pair
        
        #Nombre total d'électrons : NB=NB_UP+NB_DOWN : fixé.

        #Remplissage des états d'énergie en tenant compte du spin.
        #La somme du nombre de cases valant 1 dans tab_spin_up et dans tab_spin_down doit valoir NB.
        tab_spin_up=[0 for i in range(N)]
        tab_spin_down=[0 for i in range(N)]


        if (NB_UP%2==0):
    
            for k1 in range(N/2-NB_UP/2,N/2+NB_UP/2):
                #Electron résiduel au niveau de Fermi placé sur l'état à k<0 (convention)
                tab_spin_up[k1]=1
                
        if (NB_UP%2==1):
            
            for k1 in range(N/2-NB_UP/2,N/2+NB_UP/2+1):
                tab_spin_up[k1]=1
                
                
        if (NB_DOWN%2==0):
    
            for k1 in range(N/2-NB_DOWN/2,N/2+NB_DOWN/2):
                #Electron résiduel au niveau de Fermi placé sur l'état à k<0 (convention)
                tab_spin_down[k1]=1
                
        if (NB_DOWN%2==1):
            for k1 in range(N/2-NB_DOWN/2,N/2+NB_DOWN/2+1):
                tab_spin_down[k1]=1
 
        #print(tab_spin_down)
        #print(tab_spin_up)

        def terme_facteur_phase(n,sigma_n,m):   
            
            
            if (sigma_n==1):
                #Par convention : 1 désigne up
                res=0
                #Somme des termes de Fock sans interaction et d'Hartree sans interaction
                for j in range(N):
                    #Nombre de niveaux d'énergie parcourus : NB/2 (moitié du nombre d'électrons)
                    if (j!=n):
                        res+= -tab_spin_up[j]*cos((k_vect[n]-k_vect[j])*m*a)
                        #res+=0
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
                        res+= -tab_spin_down[j]*cos((k_vect[n]-k_vect[j])*m*a)
                        #tab_spin_down[j] vaudra 0 pour des états non occupés, 
                        #ou pour des états occupés par des spins "up" seuls ; et 1 sinon
        
                return NB-tab_spin_down[n]+res
                

        
        for n in range(N):
    
            res_spin_up=0
            res_spin_down=0
    
            for m in range(N):
                #print("terme_facteur_phase({0},{1}) = {2}".format(v,u,terme_facteur_phase(v,u)))
                res_spin_up += valeurs_I_1D[m]*terme_facteur_phase(n,1,m)
                res_spin_down += valeurs_I_1D[m]*terme_facteur_phase(n,-1,m)
        
            correction_energie_spin_up[n]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res_spin_up
            correction_energie_spin_down[n]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res_spin_down


        correction_energie=[(correction_energie_spin_up[i],correction_energie_spin_down[i]) for i in range(N)]
    
        """
        plot(k_vect,correction_energie)
        xlabel("k")
        ylabel("Delta_HF(k), polarisee en spin en eV")
        show()
        hold()
        """

        energie_corrigee_spin_up=[(E0-t0-2*t*cos(k_vect[i]*a)+correction_energie_spin_up[i]) for i in range(N)]

        comparaison_energie_up=[(energie_sans_correction[i],energie_corrigee_spin_up[i]) for i in range(N)]

        energie_corrigee_spin_down=[(E0-t0-2*t*cos(k_vect[i]*a)+correction_energie_spin_down[i]) for i in range(N)]

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
        """
    
    
        comparaison_energie_up_down=[(energie_corrigee_spin_up[i],energie_corrigee_spin_down[i]) for i in range(N)]
        
        """
        plot(k_vect,comparaison_energie_up_down) 
        #text(-10000000000, 17, r'$\Delta L_{occ}, \Delta L_{empty}$',fontsize=17)
        xlabel("k")
        ylabel("E(k)_corrigee pour UP, E(k)_corrigee pour DOWN, en eV")
        show()
        hold()
        """        
        
    
        """
        #Largeur de bandes des états occupés, largeur de bande des états vides :
        print("Largeurs de bande des états vides et occupés après correction par Hartree-Fock :")
        print("Post-correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ,Delta_vide))
        print("Post-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ/Delta_vide))
        print("Avant correction : Delta_occ= {0} eV , Delta_vide={1} eV".format(Delta_occ_1,Delta_vide_1))
        print("Avant-correction : Delta_occ/Delta_vide= {0}".format(Delta_occ_1/Delta_vide_1))
        """
    
        #print("Ces rapport doivent évoluer de façon monotone avec le remplissage NB, à N fixé")
        
        """
        #Derivée de l'énergie corrigée
        #energie_corrigee=[(E0-t0-2*t*cos(k[i]*a)+correction_energie[i]) for i in range(1000)]
        derivee_energie_corrigee=[((energie_corrigee[i+1]-energie_corrigee[i])/(k[i+1]-k[i])) for i in range((N+NB/2)/2-5,(N+NB/2)/2+5)]
        plot(k[(N+NB/2)/2-5:(N+NB/2)/2+5],derivee_energie_corrigee)
        xlabel("k")
        ylabel("dE_corr.(k)/dk en eV.m")
        #Divergence logarithmique en kF : zoomer 10 fois ne fait que doubler la taille du pic /etc...
        """
        
        
        #Calcul de l'énergie E(M) 
        res2=0 
        for k2 in range(0,N):
            res2 += tab_spin_up[k2]*energie_corrigee_spin_up[k2]
            res2 += tab_spin_down[k2]*energie_corrigee_spin_down[k2]
            #res2 += tab_spin_up[k2]*energie_sans_correction[k2]
            #res2 += tab_spin_down[k2]*energie_sans_correction[k2]
            
        energie_M[indice_aimantation]=res2
        
        
    abscisses_M=[(-NB+2*i) for i in range(0,NB+1)]
    plot(abscisses_M,energie_M)
    xlabel("Aimantation M")
    ylabel("Energie E(M) en eV")
    show()
    hold()
    
    
    
    
