# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 18:04:26 2021

Read matrix in Matrix Market Exchange Form and store it in CSR form.
If the matrix is symmetric, store it in CSC form as well.
There is also a function for the matrix to multiply with a vector.

@author: Ruei-Bo
"""


class CSR:
    def __init__(self, mat_file):
        self.FileName = mat_file
        self.shape = []
        self.val_CSR, self.col, self.row_idx, self.symmetric, self.diag = self.To_CSR(self.FileName)
        #breakpoint()
        if self.symmetric == True: 
            #breakpoint()
            self.val_CSC, self.row, self.col_idx = self.To_CSC(self.FileName)
        
    def To_CSR(self, mat_file):
        val, col, row_idx = [], [], []
        f = open(mat_file, 'r')
        has_run = False
        diag = []
        line = f.readline()
        symmetric = True if 'symmetric' in line else False
        for line in f:
            if not (line.startswith('%') or has_run):                              #filter those line begining with '%'
                shape = list(map(int, (line.split())))                               #then catch the first line (info of matrix)
                self.shape = shape
                num_sizeofmat = shape[0]
                temp = [[] for _ in range(num_sizeofmat)]                          #each row will have each list, storing all of the info of elements in the respective row
                has_run = True                                                     #just run once
            elif has_run: 
                numbers = line.split()                                             #numbers = [row_idx, col_idx, val]
                temp[int(numbers[0])-1].append(numbers)                            #
                
                
        f.close()
        for i in range(num_sizeofmat):
            for j in range(len(temp[i])):
                val.append(float(temp[i][j][2]))
                col.append(int(temp[i][j][1]))
                if j == 0:
                    row_idx.append(len(val))
                if (int(temp[i][j][0]) == int(temp[i][j][1])):
                    diag.append(float(temp[i][j][2]))
        row_idx.append(len(val)+1)
        del temp
        
        #breakpoint()
        return val, col, row_idx, symmetric, diag

    def To_CSC(self, FileName):
        val_CSC, row, col_idx = [], [], []
        f = open(FileName, 'r')
        last_col = 0
        has_run_info = False
        #diag = []
        for line in f:
            if not (line.startswith('%') or has_run_info):
                has_run_info = True
            elif has_run_info:
                numbers = line.split()
                val_CSC.append(float(numbers[2]))
                row.append(int(numbers[0]))
                if numbers[1] !=  last_col:
                    col_idx.append(len(val_CSC))
                last_col = numbers[1]
                #if numbers[0] == numbers[1]:
                #    diag.append(float(numbers[2]))
        f.close()
        col_idx.append(len(val_CSC)+1)
        
        return val_CSC, row, col_idx #, np.array(diag)
    
    def Gauss_Seidel(self, FileName):
        val, col, row_idx = [], [], []
        f = open(FileName, 'r')
        has_run = False
        diag = []
        for line in f:
            #symmetric = [True if 'symmetric' in line else False]
            if not (line.startswith('%') or has_run):                              #filter those line begining with '%'
                shape = list(map(int, (line.split())))                               #then catch the first line (info of matrix)
                self.shape = shape
                num_sizeofmat = shape[0]
                temp = [[] for _ in range(num_sizeofmat)]                          #each row will have each list, storing all of the info of elements in the respective row
                has_run = True                                                     #just run once
            elif has_run: 
                numbers = line.split()                                             #numbers = [row_idx, col_idx, val]
                temp[int(numbers[0])-1].append(numbers)                            #
       
        f.close()
        for i in range(num_sizeofmat):
            for j in range(len(temp[i])):
                if (int(temp[i][j][0]) >= int(temp[i][j][1])):
                    val.append(float(temp[i][j][2]))
                    col.append(int(temp[i][j][1]))
                    if j == 0:
                        row_idx.append(len(val))
                    if (int(temp[i][j][0]) == int(temp[i][j][1])):
                        diag.append(float(temp[i][j][2]))
        row_idx.append(len(val)+1)
        del temp
        
        return val, col, row_idx
    
    def mat_multi(self, vec):
        
        def CS_multi(val, col, row_idx, vec):
            cursor = 0
            n = len(vec)
            result = [0] * n
            for i in range(len(row_idx)-1):
                num_nonzero_ele = row_idx[i+1] - row_idx[i]
                for j in range(num_nonzero_ele):
                    result[i] += val[cursor+j]*vec[col[cursor+j]-1]
                cursor += num_nonzero_ele
            return result
        if self.symmetric == True:                                                  # if it's symmetric, (Ax + (A^T)x) - diag * x = the result of the complete matrix
            CSR_result = CS_multi(self.val_CSR, self.col, self.row_idx, vec)            # because the symmetric matrix just contains the elements of the lower triangle
            CSC_result = CS_multi(self.val_CSC, self.row, self.col_idx, vec)
            result = [a + b for a, b in zip(CSR_result, CSC_result)]
            #zero = [i for i,x in enumerate(vec) if x == 0]
            #breakpoint()
            idx = 0
            result1 = []
            for a,b in zip(result, self.diag):
                result1.append(a - b * vec[idx])
                idx += 1
            #result = [a - b for idx, a, b in enumerate(zip(result, self.diag)) if idx not in zero]
            return result1
        else:        
            result = CS_multi(self.val_CSR, self.col, self.row_idx, vec)
            #breakpoint()
            return result
        

x_nonsym = [1] * 1030
x_sym = [1] * 5357
#x_test = [0] * 3
x_test = [3,4,5]
A = CSR('aaa.mtx')
b = A.mat_multi(x_test)