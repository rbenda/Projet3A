# -*- coding: utf-8 -*-
"""
Created on Sat Oct 24 14:53:18 2015

@author: Robert
"""


get_ipython().magic(u'pylab inline')

import scipy
from scipy.integrate import quad, dblquad

import math
from math import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import numpy as np


N=100
E0=13
t0=0.5
t=2
a=math.pow(10,-10)


def integrande(z,E) :
    return 1/sqrt((1-z**2)*(1-((E0-t0-2*t*z-E)/(2*t))**2))

#Densité d'états d'énergie en 2D
def D_2D(E):
    if E >= E0-t0 and E<E0-t0+4*t:
        return N*N/(2*t*pi*pi)*quad(lambda z: integrande(z,E), -1, (E0-t0+2*t-E)/(2*t))[0] #Intégrale de -1 à alpha(E)+1
    #(alpha(E)<0)
    if E < E0-t0 and E>E0-t0-4*t:
        return N*N/(2*t*pi*pi)*quad(lambda z: integrande(z,E), (E0-t0-2*t-E)/(2*t),1)[0] #Intégrale de alpha(E)-1 à 1
    #(alpha(E)>0)
    else:
        return 0
        
        

import sympy
from sympy import Integral,limit,integrate
from sympy.abc import x
from sympy import Symbol,cos,erf,exp,sqrt


x=Symbol('x')
print(limit(1/x,x,0))

x=Symbol('x')
print(integrate(exp(-x**2)*erf(x), (x, 0, 1)))

#sympy.integrate pour les intégrales impropres
x=Symbol('x')
print(integrate(sqrt(x), (x, 0, 1)))

x=Symbol('x')
#print(integrate(D_2D(x),(x,E0-t0-4*t,E0-t0+4*t)))

