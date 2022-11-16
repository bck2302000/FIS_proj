# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 19:17:20 2021

To know how power iteration works, 
see the lec12 of FIS course slide.

@author: Bo
"""
from Read_Matrix import CSR
from Operator import *
import math
import csv
import time
import matplotlib.pyplot as plt

class PowerIte:
    
    def __init__(self, A: CSR): 
        self.A = A 
        self.z = []  
        length = self.A.shape[0]
        self.q = [[1/math.sqrt(length)] * length] # q(0), initial guess, (1/root(n)) to make L2norm be 1
        self.diff = []                            # Difference between eigenvalue and previous one
        
    
    # Solve for the largest eigenvalue.
    def powerite(self):
        k = 1           # To count the index of iteration
        Lambda = 0
        pre_Lambda = 0
        while (1):
            self.z.append(self.A.mat_multi(self.q[k-1]))           # z(k) = Aq(k-1)
            temp = vec_op(self.z[k-1], L2_norm(self.z[k-1]), '/')  # q(k) = z[k] / l2norm(z[k])
            self.q.append(temp)                                    
            Aq = self.A.mat_multi(self.q[k])                       # lambda(k) = (q(k))^H * A * q(k)
            Lambda = (multiply(self.q[k], Aq, scale = True))
            if k != 1:
                self.diff.append(abs(Lambda - pre_Lambda))         # If error < tolerance, then break  
                if self.diff[k-2] < 1e-8:
                    break
            pre_Lambda = Lambda                                    # Store the lambda as previous one
            k += 1
        self.max_Lambda = Lambda                                   # The largest eigenvalue
        
            
    # Write the data into a .txt file
    def writedata(self, elapsed_time):
        f = open('Data/power_ite_diff_'+self.A.FileName[:-4]+'.txt', 'w')  # FileName[:-4] because we don't want '.txt'
        for i in range(len(self.diff)):
            f.write(str(self.diff[i])+"\n")
        f.write("Largest EigenValue: " + str(self.max_Lambda) + '\n')
        f.write("Iterations: " + str(len(self.diff) + 1) + '\n')
        f.write("Elapsed time: " + str(elapsed_time) + 'seconds\n')        # Write the elapsed time
        f.close()
    
    def plot(self, mtxfile):
        x_axis = [i for i in range(1, len(self.diff) + 1)]                 # Set x-axis as many as the error list
        plt.plot(x_axis, self.diff)
        plt.yscale('log')
        plt.xlabel('number of iterations')
        plt.ylabel(r'$Error |\lambda^{(k)} - \lambda^{(k-1)}|$')
        plt.title('Pure Power Iteration applied on ' + mtxfile)
        plt.savefig('Data/PowerIte_' + mtxfile +'.jpg')
        plt.close()
        
    def plot_time(self, mtxfile, elapsed_time):
        x_axis = [(i*elapsed_time/len(self.diff)) for i in range(1, len(self.diff) + 1)] # Normalize it to the longest time
        plt.plot(x_axis, self.diff)
        plt.yscale('log')
        plt.xlabel('elapsed time(sec)')
        plt.ylabel(r'$Error|\lambda^{(k)} - \lambda^{(k-1)}|$')
        plt.title('Pure Power Iteration applied on ' + mtxfile)
        plt.savefig('Data/Time_PowerIte_' + mtxfile +'.jpg')
 
        
# This class is very similar to the previous one. But since tridiagonal matrix
# is not stored as CSR form, we must use different method to realize the 
# matrix-vector product. 
# The tridiagonal matrix is stored like: alpha is the diaonal, beta is the nearest 
# band beside the diagonal. 
class PowerIte_tri:
    
    def __init__(self, alpha, beta): 
        self.alpha = alpha
        self.beta = beta
        self.z = []  
        length = len(self.alpha) - 1
        self.q = [[1/math.sqrt(length)] * length]
        self.diff = []
       
    
    # Use alpha and beta to reach the matrix-vector product
    def T_multi(self, v1):
        shape = len(self.alpha) - 1
        res = []
        res.append(self.alpha[1]*v1[0] + self.beta[1] * v1[1])         # First value of the result
        # We don't put first and the last value in the for loop to prevent from 
        # accessing the elements which are out of index.
        for i in range(1,shape-1):
            temp = 0
            #e.g., beta[1]v1[0] + alpha[2]v1[1] + beta[2]v1[2]
            temp += self.beta[i] * v1[i-1] + self.alpha[i+1] * v1[i] \
                    + self.beta[i+1] * v1[i+1]
            res.append(temp)
        # Last value of the result 
        res.append(self.beta[shape-1] * v1[shape-2] +  \
                   self.alpha[shape] * v1[shape-1])
        return res
    
    
    # Do the same thing with the previous class, but different m will have different tolerance.
    def powerite(self, m):
        k = 1
        Lambda = 0
        pre_Lambda = 0
        if m < 40: 
            tol = 1e-2
        elif m >= 40 and m < 60:
            tol = 1e-4
        elif m >= 60 and m < 85:
            tol = 1e-6
        else:
            tol = 1e-10
        while (1):
            self.z.append(self.T_multi(self.q[k-1]))
            temp = vec_op(self.z[k-1], L2_norm(self.z[k-1]), '/')
            self.q.append(temp)
            Aq = self.T_multi(self.q[k])
            Lambda = (multiply(self.q[k], Aq, scale = True))
            if k != 1:
                self.diff.append(abs(Lambda - pre_Lambda))
                if self.diff[k-2] < tol:
                    break
            pre_Lambda = Lambda
            k += 1
        self.max_Lambda = Lambda
            
        
    def writedata(self, elapsed_time, mtxfile, m):
        f = open('Data/Lanczos_power_ite_diff_'+ mtxfile + '_'+ str(m)+'.txt', 'w')
        for i in range(len(self.diff)):
            f.write(str(self.diff[i])+"\n")
        f.write("Largest EigenValue: " + str(self.max_Lambda) + '\n')
        f.write("Iterations: " + str(len(self.diff) + 1) + '\n')
        f.write("Elapsed time: " + str(elapsed_time) + 'seconds\n')
        f.close()
    
    def plot(self, mtxfile, m):
        x_axis = [i for i in range(1, len(self.diff) + 1)]
        plt.plot(x_axis, self.diff)
        plt.yscale('log')
        plt.xlabel('number of iterations')
        plt.ylabel(r'$Error|\lambda^{(k)} - \lambda^{(k-1)}|$')
        plt.title('Lanczos + Power Iteration applied on ' + mtxfile +' with m = ' + str(m))
        plt.savefig('Data/Lanczos_power_' + mtxfile + '_' + str(m) +'.jpg')
        plt.close()
        
    def plot_time(self, mtxfile, m, elapsed_time):
        x_axis = [(i*elapsed_time/len(self.diff)) for i in range(1, len(self.diff) + 1)]
        plt.plot(x_axis, self.diff)
        plt.yscale('log')
        plt.xlabel('elapsed time(sec)')
        plt.ylabel(r'$Error|\lambda^{(k)} - \lambda^{(k-1)}|$')
        plt.title('Lanczos + Power Iteration applied on ' + mtxfile +' with m = ' + str(m))
        plt.savefig('Data/Time_Lanczos_power_' + mtxfile + '_' + str(m) +'.jpg')