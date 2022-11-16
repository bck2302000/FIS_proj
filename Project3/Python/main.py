# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 00:07:23 2021

@author: Bo
"""
from PowerIte import PowerIte
from PowerIte import PowerIte_tri
from Lanczos import Lanczos
from Read_Matrix import CSR
import time

# Decide execute pure power iteration or Lanzcos algorithm
def main():
    #Poweriteration()   # Pure Power Iteration
    Lanczos_Pow()       # Lanzcos Algorithm
    
# Execute pure power iteration
def Poweriteration():
    
    #mtxfile = 'nos6.mtx'                        # Smaller s.p.d. matrix
    mtxfile = 's3rmt3m3.mtx'                     # Larger s.p.d. matrix
    A = CSR(mtxfile)                             # Construct the CSR matrix form
    start = time.process_time()
    Pow = PowerIte(A)                            # Initiate the Power Iteration solver
    Pow.powerite()                               # Solve
    end = time.process_time()
    elapsed_time = end - start
    Pow.writedata(elapsed_time)                  # Write data in a .txt file
    Pow.plot(mtxfile[:-4])                       # Plot the result, error against iteration index
    Pow.plot_time(mtxfile[:-4], elapsed_time)    # Plot the result, error against time


# Execute Lanzcos algorithm
def Lanczos_Pow():
    #mtxfile = 'nos6.mtx'
    mtxfile = 's3rmt3m3.mtx'
    A = CSR(mtxfile)
    m = 50                                       # Dimension of Krylov space. Users can adjust it
    start = time.process_time()
    Lan = Lanczos(A, m)                          # Initialize Lanczos solver
    Lan.T_construct()                            # Construct the tridiagonal matrix
    Lan_Pow = PowerIte_tri(Lan.alpha, Lan.beta)  # Use power iteation to solve the tridiagonal matrix
    Lan_Pow.powerite(m)                          # Solve for the largest eigenvalue
    end = time.process_time()
    elapsed_time = end - start
    Lan_Pow.writedata(elapsed_time, mtxfile[:-4], m)  # Write data in a .txt file
    Lan_Pow.plot(mtxfile[:-4], m)                     # Plot the result, error against iteration index
    Lan_Pow.plot_time(mtxfile[:-4], m, elapsed_time)  # Plot the result, error against time
    
if __name__ == "__main__":
    main()