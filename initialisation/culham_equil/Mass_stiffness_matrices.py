# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 08:48:06 2015

@author: GF245157
"""

from math import *
import scipy as sc
import sympy as sp
import numpy as np
from fractions import *

def mass_matrix(u0,v0,n):
    """
        Computation of \chi_{n,n,n} (u0 + n, v0 + n)
    """
    u = np.copy(u0)
    v = np.copy(v0)
    x = u + v
    sqrt3y = v - u
    x = - abs(x)
    sqrt3y = abs(sqrt3y)
    u = (x - sqrt3y)/2
    v = (x + sqrt3y)/2
    
    id=np.nonzero(v>0)
    v[id]=-v[id]
    u[id]=u[id]+v[id]

    id=np.nonzero(v>u/2)
    v[id]=u[id]-v[id]
#    print(u)
#    print(v)
    u = u.astype(int)
    v = v.astype(int)
#    print(u)
#    print(v)

    val=np.zeros(u.shape, dtype=Fraction)
    print(val)

    for K in range(-n, int(ceil(np.max(u)))):
        K0=Fraction(K)
        for L in range(-n, int(ceil(np.max(v)))):
            L0=Fraction(L)
            for i in range(0, min(n+K, n+L)+1):
                coeff=Fraction((-1)**(K+L+i)*binomialCoefficient(n,i-K)*binomialCoefficient(n,i-L)*binomialCoefficient(n,i))
#                print('coeff = ' + str(coeff))
                for d in range(0,n):
                    aux=abs(v-L0-u+K0)
#                    print('aux = ')
#                    print(aux)
                    aux2=(u-K0+v-L0-aux)/2
                    aux2[np.nonzero(aux2<0)]=0
#                    print('aux2 = ')
#                    print(aux2)
                    val=val+coeff*Fraction(binomialCoefficient(n-1+d,d),factorial(2*n-1+d)*factorial(n-1-d))*aux**(n-1-d)*aux2**(2*n-1+d)
                #end for
            #end for
        #end for
    #end for
                    
    return val
#end def

def stiffness_matrix(u0,v0,n):
    """
        Computation of
        (\partial_{r_1, r_1} + \partial_{r_2, r_2} - \partial_{r_1, r_2}) \chi_{n,n,n} (u_0 + n, v_0 + n)
    """
    u = np.copy(u0)
    v = np.copy(v0)
    x = u + v
    sqrt3y = v - u
    x = - abs(x)
    sqrt3y = abs(sqrt3y)
    u = (x - sqrt3y)/2
    v = (x + sqrt3y)/2
    
    id=np.nonzero(v>0)
    v[id]=-v[id]
    u[id]=u[id]+v[id]

    id=np.nonzero(v>u/2)
    v[id]=u[id]-v[id]
#    print(u)
#    print(v)
    u = u.astype(int)
    v = v.astype(int)
#    print(u)
#    print(v)

    val=np.zeros(u.shape, dtype=Fraction)

    for K in range(-n, int(ceil(np.max(u)))):
        K0=Fraction(K)
        for L in range(-n, int(ceil(np.max(v)))):
            L0=Fraction(L)
            for i in range(0, min(n+K, n+L)+1):
                coeff=Fraction((-1)**(K+L+i))*binomialCoefficient(n,i-K)*binomialCoefficient(n,i-L)*binomialCoefficient(n,i)
                for d in range(0, n):
                    aux=abs(v-L0-u+K0)
                    aux2=(u-K0+v-L0-aux)/2
                    aux2[np.nonzero(aux2<0)]=Fraction(0)
                    val=val+Fraction(binomialCoefficient(n-1+d,d),factorial(2*n-3+d)*factorial(n-1-d))*coeff*aux**(n-1-d)*aux2**(2*n-3+d)
                    if n-1-d > 0:
                        val=val-Fraction(3)*coeff*Fraction(binomialCoefficient(n-1+d,d),factorial(2*n-2+d)*factorial(n-2-d))*aux**(n-2-d)*aux2**(2*n-2+d)
                    #end if
                    if n-2-d > 0:
                        val=val+Fraction(3)*coeff*Fraction(binomialCoefficient(n-1+d,d), factorial(2*n-1+d)*factorial(n-3-d))*aux**(n-3-d)*aux2**(2*n-1+d)
                    #end if
                #end for
            #end for
        #end for
    #end for
    
    return val
#end def
    
def binomialCoefficient(n, k):
    """
        Computation of the binomial coefficient (n, k)
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    k = min(k, n - k) # take advantage of symmetry
    c = Fraction(1)
    for i in range(k):
        c = c * (n - i) / (i + 1)
    return int(round(c))
#end def

def mass_stiff_matrices(n0):
    """
        Computation of the coefficients which will be found in the mass
        and stiffness matrices. The result is two matrices (let's call them MM
        and SM)
        
        For a BS centered in (k_{1,i},k_{2,i}) and another one centered in
        (k_{1,j},k_{2,j}), the corresponding coefficient for the mass matrix
        (resp. for the stiffness matrix) at position (i,j) is the coefficient
        of MM (resp. SM) at position:
        (2*n0 + (k_{1,i} - k_{1,j})/h, 2*n0 + (k_{2,i} - k_{2,j})/h)
    """
    n = 2*n0
    u = range(-n+1, n)
    v = range(-n+1, n)
    U, V = np.meshgrid(u,v)
    U = np.array(np.round(U), dtype=Fraction)
    V = np.array(np.round(V), dtype=Fraction)
    print(U)
    print(V)
    
    mM = mass_matrix(U,V,n)
    sM = stiffness_matrix(U,V,n)
    
    return sp.Matrix(mM), sp.Matrix(- sM)
