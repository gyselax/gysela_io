# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:18:39 2015

@author: GF245157
"""

from numpy import zeros, shape
from scipy.sparse import csr_matrix

def coefbinbiz(N):
    """
        Compute the coefficients found in the h-refinement relation
        with their induction relation
    """
    n=N+1
    A=zeros((2*N+5,2*N+5,n), dtype=int)

    A[2,2,0] = 1

    for n0 in range(1,n):
        for i in range(2, (2*n0+3)):
            for j in range(2,(2*n0+3)):
                A[i,j,n0]=2*A[i-1,j-1,n0-1]+A[i,j,n0-1]+A[i,j-1,n0-1]+A[i-1,j,n0-1]+A[i-2,j-1,n0-1]+A[i-1,j-2,n0-1]+A[i-2,j-2,n0-1]
            #end for
        #end for
    #end for

    return A[2:(2*N+3),2:2*N+3,n-1]

def gth_to_htg(gth_BS):
    """
        From the np array which gives the hexagonal coordinates
        of the center with the position of the BS, return a dictionary which
        gives the global coordinate i of the BS when we give it the center of the BS
    """
    htg = {}
    for i in range(0, gth_BS.shape[0]):
        htg[gth_BS(i, 0), gth_BS(i, 1)] = i
    #end for
    
    return htg
#end def


def h_refinement(htg_BS_h, gth_BS_2h, deg_BS, h):
    """
        From the array 'global_to_hex' of the mesh with step 2*h (the center
        of the Box-splines of the 2*h-step mesh at index i must be given by
        (gth_BS_2h[i,0], gth_BS_2h[i,1]) in hexagonal coordinates) and
        the dictionary 'hex_to_global' of the mesh with step h (the index of
        the Box-splines of the h-step mesh which center is (k1, k2)
        in hexagonal coordinates must be given by htg_BS_h[k1,k2] ), the degree
        of the Box-splines and h, return the h-refinement matrix.
        
        This works only for Box-splines of type I (hexagonal or not):
        the Box-splines must have 3 directions u_1, u_2 and u_3 = u_1 + u_2
        and the same degree for every direction
    """
    M = csr_matrix((len(htg_BS_h), gth_BS_2h.shape[0]), dtype=int)
    coeff = coefbinbiz(deg_BS)
    
    for i in range(0, gth_BS_2h.shape[0]):
        k1 = gth_BS_2h[i, 0]
        k2 = gth_BS_2h[i, 1]
        for lambd in range(0,2*deg_BS+1):
            for gamma in range(0,2*deg_BS+1):
                if coeff[gamma, lambd] != 0:
                    if (2*k1+gamma-deg_BS, 2*k2+lambd-deg_BS) in htg_BS_h.keys():
                        M[htg_BS_h[2*k1+gamma-deg_BS, 2*k2+lambd-deg_BS], i] = coeff[gamma, lambd]
                    #end if
                #end if
            #end for
        #end for
    #end for
    
    return 2**(2-3*deg_BS)*M
