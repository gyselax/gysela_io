# -*- coding: utf-8 -*-
"""
Created on Thu May 28 14:33:38 2015

@author: GF245157
"""

from math import *
import sympy as sp
import numpy as np
from fractions import *
import matplotlib.pyplot as mpp
import cmath

def ringradiussquare(n_max):
    """
        Return 2 results
        radius is a dict which gives the 2-norm to the square of vectors n1 r_1 + n2 r_2 with n1 and n2 integers,
        such that M = n1+n2 and N = n1-n2. This value is then given by
        (N**2 + 3 * M**2)/4. We take points such that this value is less than n_max.
        Rradius is a sorted list which gives all the radius which are in radius
    """
    interm=set()
    radius={}
    for N in range(-2*n_max, 2*n_max+1):
        a=int(math.sqrt((4*n_max**2-N**2)/3))
        mod = N % 2
        for x in range(-int((a+mod)/2), int((a-mod)/2)+1):
            M = mod + 2*x
            r = round((N**2+3*M**2)/4)
            interm = interm | set([r])
            radius[N,M] = r
        #end for
    #end for
    Rradius=sorted(list(interm))
        
    return Rradius, radius
#end def

def ringnumber(Rradius, radius, K):
    """
        Return the circle number of all points of radius s.t. its circle number is less than K
    """
    Rnumber={}
    for ((N,M),r) in radius.items():
        numb = Rradius.index(r)
        if numb <= K:
            Rnumber[N,M] = numb
        #end if
    #end for
            
    return Rnumber
#end def

def make_coeff_function(Rnumber, coeff, n):
    """
        Return the coefficients of the Taylor expansion of p(\omega_1, 0)
        depending on h[0], ..., h[n] given in coeff.
        coeff_function[i] corresponds to the coefficient of omega**i
    """
    coeff_function = sp.zeros(1, n+1)
    for (N,M) in Rnumber.keys():
        if (N,M) == (0,0):
            coeff_function[0] = coeff_function[0] + coeff[0]
        elif (N,M) > (0,0):
            # Each point is linked with its opposite in order to have a cosinus
            for i in range(0, n+1):
                coeff_function[i] = coeff_function[i] + 2*(-1)**i*(Fraction((N+M),2))**(2*i)/factorial(2*i)*coeff[Rnumber[N,M]]
            #end for
        #end if
    #end for
    
    return coeff_function
#end def

def make_function(Rnumber, coeff_function_modif, n_modif):
    """
        Build the function of the Taylor expansion of p(\omega_1, 0)
        with the coefficients coeff_function_modif until the order n_modif
    """
    global omega
    omega = sp.Symbol('omega')
    p = 0
    for i in range(0,n_modif+1):
        p = p + coeff_function_modif[i]*omega**(2*i)
    #end for
#    print(p)
    return p
#end def

def makelambd_2(n):
    """
        Compute the coefficients of the Taylor expansion of B_n(\omega_1, 0) until the order n
    """
    lambd=sp.zeros((n+1,1))
    j=np.zeros(n, dtype=int, order='C')
    while j[0]!=n:
        increase(j, n)
        s=sum(j)
        if s<=n:
            coeff=Fraction(1)
            for i in range(0,n):
                coeff=coeff/Fraction(factorial(2*j[i]+1))
            #end for
            lambd[s]=lambd[s]+coeff
        #end if
    #end while
            
    for i in range(1,n+1):
        lambd[i]=lambd[i]*Fraction(1, 4**i)*(-1)**i
    #end for
        
    lambd[0]=Fraction(1)
    
    return lambd
#end def

def increase(j, n):
    """
        Increase the counters in makelambd_2
    """
    if sum(j) == n:
        i=-1
        while j[i] == 0:
            i=i-1
        #end while
        j[i]=0
        j[i-1]=j[i-1]+1
    else:
        j[-1]=j[-1]+1
#end def

def coeff_computation_2(n):
    """
        Principle function which builds the functions, the equations and solve them
    """
    Rradius, radius = ringradiussquare(n)
    Rnumber = ringnumber(Rradius, radius, n)
    lambd = makelambd_2(2*n)
    
    # Plot of the taken points
    x = [N[1]/2 for N in Rnumber.keys()]
    y = [N[0]*sqrt(3)/2 for N in Rnumber.keys()]
    mpp.plot(x, y, 'o')
    
    # Definition of the symbolic coefficients h[i], in coeff
    coeff = {}
    for i in range(0,n+1):
        coeff[i] = sp.Symbol("".join([' h', str(i), ' ']))
    #end for
    
    # Build of the coefficients of the Taylor expansion of p(\omega_1, 0)
    coeff_function = make_coeff_function(Rnumber, coeff, n)
    coeff_function_modif = coeff_function[:,:]
    
    # Matrix and vector of the linear equations: not needed if we use "solve"
#    M = sp.zeros((n+1, n+1))
#    V = sp.zeros((n+1, 1))
    
    # Matrix which countains the linear equations
    Lineq = sp.zeros(1, n+1)
    
    # Beginning of the loop which builds the equations
    for i in range(0,n+1):
        # Computation of the i^th coefficient in the Taylor expansion of 1/p(\omega_1, 0)
        # taking into account of the previous equations
        p = make_function(Rnumber, coeff_function_modif, i)
        b = 1/p
        result=sp.diff(b, omega, 2*i).subs(omega, 0)/factorial(2*i)
        
        # for the first iteration, the equation we have is 1/(linear equation in h[0], ..., h[n]) = 1
        if i == 0:
            result = 1/result
        #end if
        
        # Computation of the coefficients of M and V
#        gamma = result
#        for k in range(0, n+1):
#            gamma = gamma.subs(coeff[k], 0)
#        #end for
#        V[i] = lambd[i] - gamma
#        for k in range(0, n+1):
#            result_interm = result
#            for j in range(0, n+1):
#                if j == k:
#                    result_interm = result_interm.subs(coeff[k],1)
#                else:
#                    result_interm = result_interm.subs(coeff[j],0)
#                #end if
#            #end for
#            result_interm = result_interm - gamma
#            M[i,k] = result_interm
#        #end for
#        
#### ####        result_interm = result.subs(coeff[k], 1).subs(coeff.keys(), np.zeros(n+1))
        
        # Linear equation put into the matrix
        Lineq[i] = result-lambd[i]
        
        # Modification of the coefficients by taking into account the new equation
        coeff_function_modif[i] = coeff_function_modif[i].subs(coeff[1], sp.solve(result-lambd[i], coeff[1])[0])
    #end for
    
    # Computation of the results with matrix and vector
#    print(M)
#    results = M.inv()*V
    
    # Computation of the results with "solve"
    results_1 = sp.solve(Lineq, list(coeff.values()), dict=True)
    print(' Results : ', results_1)
    if results_1 == []:
        print('No solution')
    else:
        results = sp.zeros(1,n+1)
        for i in range(0,n+1):
            results[i] = results_1[0][coeff[i]]
        #end for
    
        coeff_validation(n, Rnumber, results)
    
        return results
    #end if
#end def
    
def coeff_validation(n, Rnumber, results):
    """
        Print the Taylor expansion of B_n(\omega_1,0) - 1/p(\omega_1,0)
        at the order 2*(n+1) with the computed coefficients to valide them
    """
    B = (sp.sin(omega/2)*2/omega)**(2*n)
    b = 0
    for (N,M) in Rnumber.keys():
        if (M,N) > (0,0):
            b = b + 2 * results[Rnumber[N,M]] * sp.cos(((M+N)/2)*omega)
        elif (M, N) == (0, 0):
            b = b + results[0]
        #end if
    #end for
    print((B-1/b).series(omega, 0, 2*(n+1)))
#end def