# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 18:04:26 2021

@author: Bo
"""

import numpy as np

def To_CSR(mat_file):
    val, col, row_idx = [], [], []
    f = open(mat_file, 'r')
    has_run = False
    symmetric = True
    for line in f:
        if not (line.startswith('%') or has_run):                              #filter those line begining with '%'
            num_sizeofmat = int(line.split()[0])                               #then catch the first line (info of matrix)
            temp = [[] for _ in range(num_sizeofmat)]                          #each row will have each list, storing all of the info of elements in the respective row
            has_run = True                                                     #just run once
        elif has_run: 
            numbers = line.split()                                             #numbers = [row_idx, col_idx, val]
            temp[int(numbers[0])-1].append(numbers)                            #
            if int(numbers[0]) < int(numbers[1]):                              #if any row_idx < col_idx, this matrix is symmetry
                symmetric = False
            
    f.close()
    for i in range(num_sizeofmat):
        for j in range(len(temp[i])):
            val.append(float(temp[i][j][2]))
            col.append(int(temp[i][j][1]))
            if j == 0:
                row_idx.append(len(val))
    row_idx.append(len(val)+1)
    del temp
    
    return val, col, row_idx, symmetric

def To_CSC(mat_file):
    val, row, col_idx = [], [], []
    f = open(mat_file, 'r')
    last_col = 0
    has_run_info = False
    diag = []
    for line in f:
        if not (line.startswith('%') or has_run_info):
            has_run_info = True
        elif has_run_info:
            numbers = line.split()
            val.append(float(numbers[2]))
            row.append(int(numbers[0]))
            if numbers[1] !=  last_col:
                col_idx.append(len(val))
            last_col = numbers[1]
            if numbers[0] == numbers[1]:
                diag.append(float(numbers[2]))
    f.close()
    col_idx.append(len(val)+1)
    
    return val, row, col_idx, np.array(diag)

def mat_multi(mat_file, vec):
    
    def CS_multi(val, col, row_idx, vec):
        cursor = 0
        result = np.zeros_like(vec)
        for i in range(len(row_idx)-1):
            num_nonzero_ele = row_idx[i+1] - row_idx[i]
            for j in range(num_nonzero_ele):
                result[i] += val[cursor+j]*vec[col[cursor+j]-1]
            cursor += num_nonzero_ele
        return result
    
    val_CSR, col, row_idx, symmetric = To_CSR(mat_file)
    if symmetric:                                                              # if it's symmetric, (Ax + (A^T)x) - diag = the result of the complete matrix
        val_CSC, row, col_idx, diag = To_CSC(mat_file)                         # because the symmetric matrix just contains the elements of the lower triangle
        CSR_result = CS_multi(val_CSR, col, row_idx, vec)
        CSC_result = CS_multi(val_CSC, row, col_idx, vec)
        result = CSR_result + CSC_result
        result = result + diag
        return result
    else:
        result = CSR_result = CS_multi(val_CSR, col, row_idx, vec)
        return result
    
    
def Gram_Schmidt(A, r0, m):
    

x_nonsym = np.ones(1030)
x_sym = np.ones(5357)
#b = mat_multi('orsirr_1.mtx', x_nonsym)
b = mat_multi('s3rmt3m3.mtx', x_sym)