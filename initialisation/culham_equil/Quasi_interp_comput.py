# -*- coding: utf-8 -*-
"""
Created on Tue May 26 15:55:09 2015

@author: GF245157
"""

import numpy as np
import fractions as fr
import math
import sympy as sp

def makematrix(n):
    M=sp.zeros((n+1,n+1))
    for i in range(0,n+1):
        M[0,i]=fr.Fraction(2)
        for j in range(1,n+1):
            M[j,i]=(-1)**j*fr.Fraction(2*i**(2*j),math.factorial(2*j))
        #end for
    #end for
    M[0,0]=fr.Fraction(1)
    
    return M
#end def

def makelambd(n):
    lambd=sp.zeros((n+1,1))
    j=np.zeros(n, dtype=int, order='C')
    while j[0]!=n:
        increase(j, n)
        s=sum(j)
        if s<=n:
            coeff=fr.Fraction(1)
            for i in range(0,n):
                coeff=coeff/fr.Fraction(math.factorial(2*j[i]+1))
            #end for
            lambd[s]=lambd[s]+coeff
        #end if
    #end while
            
    for i in range(1,n+1):
        lambd[i]=lambd[i]*fr.Fraction(1, 4**i)*(-1)**i
    #end for
        
    lambd[0]=fr.Fraction(1)
    
    return lambd
#end def

def increase(j, n):
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

def coeff_computation(n):
    M=makematrix(n)
#    V=makevector(n)
    V=makelambd(n)
    
    return M.inv()*V
#end def