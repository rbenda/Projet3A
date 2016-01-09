# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 18:12:31 2015

@author: Robert
"""

"""
Ce code permet de calculer la correction à l'énergie d'un niveau électronique donné (pour un réseau 2D),
en écrivant les intégrales I_2D(l1,p1) en coordonnées sphériques et en tirant les variables radiales
selon des lois gaussiennes et les variables angulaires uniformément.
L'invariance par translation dans la somme a été prise en compte.
"""

get_ipython().magic(u'pylab inline')
import numpy 
import scipy
import math
from math import * 
import random
import time


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


e2=2.3*math.pow(10,-28) 
N=30
E0=13
t0=0.5
t=2
a=math.pow(10,-10)

l2=0
p2=0

#d : inférieur ou égal à a/4. pour une gaussienne (ainsi 1 % de recouvrement des gaussiennes
# de 2 orbitales atomiques localisée voisines en a/2)
d=a/4.
#Tester d'autres d ??





#Nombre d'ELECTRONS mis dans le système
#On remplit les états de façon non magnétique : en partant de l'état fondamental ; 
#2 électrons (spins up et down) par état

#Si le nombre d'électrons est impair, on a forcément un état d'énergie avec un électron
#uniquement et donc un seul spin représenté. La correction calculée sera polarisée en spin.

#Nombre de valeurs de k_x occupées 
alpha1=30
#Nombre de valeurs de k_y occupées 
alpha2=30

NB=alpha1*alpha2
#NB maximal : 2*(N^2)

#Demi-remplissage : NB=N 
#(N niveaux en tout ; 2 électrons par niveaux : donc 2N électrons au total au maximum)

k_x=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]

k_y=[(-(2*pi/(N*a))*floor(N/2.)+i*(2*pi/(N*a))) for i in range(0,N)]


def distance_reseau_2D(x,y,z,x1,y1,z1):
    #§Minimum entre 4 éléments
    var1=sqrt((x-(x1-N*a))**2+(y-y1)**2+(z-z1)**2)
    var2=sqrt((x-x1)**2+(y-(y1-N*a))**2+(z-z1)**2)
    var3=sqrt((x-(x1-N*a))**2+(y-(y1-N*a))**2+(z-z1)**2)
    return min(var1,var2,var3,sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)) 


def terme_facteur_phase(n_x,n_y,l1,p1):
    res=0

    for j_x in range(alpha1):
        for j_y in range(alpha2):
            if ((N-alpha1)/2+j_x!=n_x) or ((N-alpha2)/2+j_y!=n_y):
                #rajouter la condition : j est un état occupé de même spin ? 
                #Faire une fonction qui dit si l'état est occupé ou non ?
                res+= -cos((k_x[n_x]-k_x[(N-alpha1)/2+j_x])*(l1-l2)*a+(k_y[n_y]-k_y[(N-alpha2)/2+j_y])*(p1-p2)*a)
    
    #if (n_x >= (N-alpha1)/2) and (n_x < (N+alpha1)/2) and (n_y >= (N-alpha2)/2) and (n_y < (N+alpha2)/2):
    #    return NB+res
    #else 
        #(k_n_x,k_n_y) n'est pas un état occupé
    #    return NB-1+res
    
    return res
    #Effet du terme de Fock seul

    
def F(l1,p1,r1,r2,theta1,theta2,phi1,phi2):
    #return (pi**2)*(2*pi)**2*(pi*(d/a)**2/2)*r1**2*r2**2*sin(theta1)*sin(theta2)/sqrt(r1**2+r2**2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))+2*l1*(r2*sin(theta2)*cos(phi2)-r1*cos(phi1)*sin(theta1))+l1**2+2*p1*(r1*sin(theta1)*sin(phi1)-r2*sin(theta2)*sin(phi2))+p1**2)
    return r1**2*r2**2*sin(theta1)*sin(theta2)/distance_reseau_2D(r1*sin(theta1)*cos(phi1)-l1,r1*sin(theta1)*sin(phi1)-p1,r1*cos(phi1),r2*sin(theta2)*cos(phi2),r2*sin(theta2)*sin(phi2),r2*cos(phi2))
#Facteur exp(-2*r1**2/(d/a)**2)*exp(-2*r1**2/(d/a)**2) dans l'intégrande : contenu dans
#le tirage de r1 et r2 selon des gaussiennes    
#Premier facteur : facteur de renormalisation de la densité de la gaussienne
    
nb=1000.

"""
r1_0=0.
r2_0=0.
theta1_0=0.
theta2_0=0.
phi1_0=0.
phi2_0=0.
"""
 
def I_2D(l1,p1):
    
    Y=[0 for i in range(6)]    
    res=0

    for i in range(999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
    
        
        Y[0]=numpy.random.normal(0,(d/a)/2.)
        Y[1]=numpy.random.normal(0,(d/a)/2.)
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(0,2*pi)
        Y[5]=random.uniform(0,2*pi)
        
        res+= F(l1,p1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
   
    return (pi**2/d)*(a/d)**3*(res/nb)
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


"""
#Avec cette méthode, on recalcule pour chaque m 
#I_1D(m) avec nb tirages du 6-uplet (rho1,rho2,theta1,theta2,phi1,phi2)
def Delta_Hartree_Fock(n_x,n_y):
    res2=0
    tab_I_2D=[[0 for i in range(N)] for j in range(N)]
    
    for l1 in range(N):
        for p1 in range(N):
            tab_I_2D[l1][p1]=I_2D(l1,p1)*terme_facteur_phase(n_x,n_y,l1,p1)
            res2+= tab_I_2D[l1][p1]
        #Le calcul de I_1D(m) nécessite à chaque m nb tirages...
        
    #abscisses_m=[i for i in range(N)]
    #plot(abscisses_m,tab_I_1D)
    return (e2/N)*res2
    
#print(Delta_Hartree_Fock(2))


nb_tirages=10000.

#Avec cette méthode, on évalue pour chaque tirage la valeur de l'intégrande de I_1D
#pour chaque m.
#Rajouter variable de spin en argument ?
def Delta_Hartree_Fock_bis(n_x,n_y):
    res=[[0 for i in range(N)] for j in range(N)]
    
    Y=[0 for i in range(6)]

    for i in range(9999):
        #Tirage d'un (i+1)ème uplet de variables aléatoires
        Y[0]=numpy.random.normal(0,(d/a)/sqrt(2.))
        Y[1]=numpy.random.normal(0,(d/a)/sqrt(2.))
        Y[2]=random.uniform(0,pi)
        Y[3]=random.uniform(0,pi)
        Y[4]=random.uniform(-pi,pi) 
        Y[5]=random.uniform(-pi,pi)
     
        #tab=[0 for i in range(N)]
        
        for l1 in range(N):
            for p1 in range(N):
                res[l1][p1] += F(l1,p1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
                #tab[m] = F(m,p1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5])
            
        #Vérification que à p1 fixé, pour un tirage des variables donné, l1 -> F(l1,p1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]) évolue bien en 1/|l1|
        #abscisses_m=[i for i in range(10)]
        #print(Y)
        #print(F(1,1,Y[0],Y[1],Y[2],Y[3],Y[4],Y[5]))
        #plot(abscisses_m,tab[0:10])
        
    I_2D_estime=[[0 for i in range(N)] for j in range(N)]
    
    for l1 in range(N):
            for p1 in range(N):
                I_2D_estime[l1][p1]=res[l1][p1]/nb_tirages
    
    resultat_correction_energie=0
    for l_1 in range(N):
        for p_1 in range(N):
            resultat_correction_energie += I_2D_estime[l1][p1]*terme_facteur_phase(n_x,n_y,l_1,p_1)
        
    return (e2/N)*resultat_correction_energie
"""



"""   
n=5
 
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

"""
tab=[[0 for i in range(N)] for j in range(N)]
for m in range(N):
    for q in range(N):
        tab[m][q]=Delta_Hartree_Fock_bis(m,q)
        print("Delta_Fock_w_s_i_(k_x[{0}]={1},k_y[{2}]={3}) : {4} ".format(m,k_x[m],q,k_y[q],tab[m][q]))
"""

#numpy.savetxt('Correction_energie_2D_bis.txt',tab,newline='\n')


valeurs_I_2D=[[0 for i in range(N)] for j in range(N)]
#On calcule les valeurs de I_1D une fois pour toutes (en faisant beaucoup de tirages pour un m donné et ainsi avoir peu d'incertitudes)
#(elles sont indépendantes de n ; n n'intervient que dans le facteur de phase)
for l1 in range(N):
    for p1 in range(N):
        valeurs_I_2D[l1][p1]=I_2D(l1,p1)
    

#print(valeurs_I_2D)

null=linspace(0.,0.,N)
correction_energie=np.outer(null,np.ones(np.size(null)))
    
for n_x in range(N):
    for n_y in range(N):
        #tps1=time.clock()
    
        res=0
        for l1 in range(N):
            for p1 in range(N):
                #print("terme_facteur_phase({0},{1}) = {2}".format(v,u,terme_facteur_phase(v,u)))
                res+= valeurs_I_2D[l1][p1]*terme_facteur_phase(n_x,n_y,l1,p1)
        
        correction_energie[n_x][n_y]=(1/(1.6*math.pow(10,-19)))*(e2/N)*res
    
        #☺tps2=time.clock()
        #print("Delta_Hartree_Fock({0},{1})= {2}  Temps calcul : {3}".format(n_x,n_y,correction_energie[n_x][n_y],tps2-tps1))


#Graphe en 3D de la correction de l'énergie par le terme de Hartree Fock
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D


#plt.show() 

k_x_1=np.outer(k_x,np.ones(np.size(k_x)))

k_y_1=np.outer(np.ones(np.size(k_y)),k_y)


energie_2D_corrigee=[[(E0-t0-2*t*(cos(k_x[i]*a)+cos(k_y[j]*a))+correction_energie[i][j]) for i in range(N)] for j in range(N)]
#energie_corrigee=[[(energie_2D[i][j],correction_energie[i][j]) for i in range(N)] for j in range(N)] 
#energie_corrigee=[[(energie_2D[i][j]+correction_energie[i][j]) for i in range(N)] for j in range(N)]

ax.plot_surface(k_x_1, k_y_1,correction_energie,rstride=10, cstride=10, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
#La fonction plot_surface trace à partir de 3 matrices A, B et C, l'ensemble des points 
#de coordonnées (A[i][j], B[i][j], C[i][j]) et les relie pour former une surface.

plt.title("Correction given by the Fock term (without self-interaction) for the 2D lattice,  in eV", fontsize=10)
plt.show()
plt.hold()
#ax.contour(k_x_1, k_y_1, energie_corrigee,zdir='z')
 

#ax.set_zlim(E0-t0-4*t, E0-t0+4*t)

ax.zaxis.set_major_locator(plt.LinearLocator(10))
ax.zaxis.set_major_formatter(plt.FormatStrFormatter('%.02f'))

X,Y= meshgrid(k_x,k_y)
pcolor(X,Y,correction_energie)
show()


#fig.colorbar(surf, shrink=0.5, aspect=10)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d') #Créaton d'axes 3D

ax.plot_surface(k_x_1, k_y_1,energie_2D_corrigee,rstride=10, cstride=10, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
plt.show()
plt.hold()

"""
X,Y= meshgrid(k_x,k_y)
pcolor(X,Y,energie_2D)
show()
"""

