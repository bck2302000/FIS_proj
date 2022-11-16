# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 09:05:34 2021

To get better idea about Lanzcos algorithm,
see the lecture12 slide of FIS course.

@author: Bo
"""
from Read_Matrix import CSR
from Operator import *
import math

class Lanczos:
    
    def __init__(self, A, m):
        self.A = A
        length = self.A.shape[0]
        self.v = [[0] * length, [1/math.sqrt(length)] * length]
        self.beta = [0]
        self.alpha = [0]
        self.m = m     # Dimension of Krylov space
    
    
    # To construct the tridiagonal matrix
    def T_construct(self):
        for j in range(1,self.m+1):
            
            # w = Av[j] - beta[j-1]v[j-1]
            w = vec_subtract(self.A.mat_multi(self.v[j]), \
                             vec_op(self.v[j-1], self.beta[j-1], '*'))
                
            self.alpha.append(multiply(self.v[j], w, scale = True))     # alpha[j] = (v^T)w
            w = vec_subtract(w, vec_op(self.v[j], self.alpha[j], "*"))  # w = w - alpha[j]v[j]
            self.beta.append(L2_norm(w))                                # beta[j] = L2norm(w)
            self.v.append(vec_op(w, self.beta[j], '/'))                 # v[j+1] = w / beta[j]
    
        