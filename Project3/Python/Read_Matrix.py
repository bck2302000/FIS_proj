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
        self.val_CSR, self.col, self.row_idx, self.symmetric, self.diag \
            = self.To_CSR(self.FileName)  # Convert the matrix from the matrix market file to CSR form.
        
        # If the matrix market file specifies it is a symmetric matrix, 
        # convert to CSC form as well. The matrix-vector product for 
        # symmetric matrix can be realized by A^T * v + A * v - diag.
        if self.symmetric == True: 
            self.val_CSC, self.row, self.col_idx = self.To_CSC(self.FileName)
        
    
    # Convert COO to CSR form
    def To_CSR(self, mat_file):
        val, col, row_idx = [], [], []
        f = open(mat_file, 'r')        # Open the matrix file as f
        has_run = False                # Whether the matrix info line has been read?
        diag = []                      # Store the value on the diagonal
        line = f.readline()
        symmetric = True if 'symmetric' in line else False   # if it's a symmetric matrix, remark it
        
        # Read the matrix market file line by line. Put them in the list first.
        for line in f:
            if not (line.startswith('%') or has_run):      # Filter those line begining with '%'
                shape = list(map(int, (line.split())))     # then catch the first line (info of matrix)
                self.shape = shape
                num_sizeofmat = shape[0]
                temp = [[] for _ in range(num_sizeofmat)]  # Each row has each list, storing all of the info of elements in the respective row
                has_run = True                             # We just read the matrix info once.
            elif has_run: 
                numbers = line.split()                     # Numbers = [row_idx, col_idx, val]
                temp[int(numbers[0])-1].append(numbers)    # Store the info into the respective list.
        f.close()                                          # Close the file
        
        
        # Based on the temp list, which is ordered by the row, we can store them
        # in CSR form
        for i in range(num_sizeofmat):           # Size of the matrix
            for j in range(len(temp[i])):        # The number of non-zero val in each row
                val.append(float(temp[i][j][2])) # Val and column can be directly appended
                col.append(int(temp[i][j][1]))
                if j == 0:                       # When meeting the first non-zero ele in the row,
                    row_idx.append(len(val))     # append row_idx. For detailed please google CSR form
                if (int(temp[i][j][0]) == int(temp[i][j][1])):   # Store the diagonal
                    diag.append(float(temp[i][j][2]))
        row_idx.append(len(val)+1)               # Append the (number of val + 1) in the last row_idx 
        del temp
        return val, col, row_idx, symmetric, diag


    # Convert COO to CSC form. This is to tackle the matrix-vector product of symmetric matrix
    # The method is similar to convert to CSR, but some order are different
    def To_CSC(self, FileName):
        val_CSC, row, col_idx = [], [], []
        f = open(FileName, 'r')
        last_col = 0
        has_run_info = False
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
        f.close()
        col_idx.append(len(val_CSC)+1)
        
        return val_CSC, row, col_idx
    
    def Gauss_Seidel(self, FileName):
        val, col, row_idx = [], [], []
        f = open(FileName, 'r')
        has_run = False
        diag = []
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
    
    
    # Perform matrix-vector product.
    def mat_multi(self, vec):
        
        def CS_multi(val, col, row_idx, vec):
            cursor = 0                            # To record which values we are computing
            n = len(vec)                          # The size of the vecotr
            result = [0] * n                      # Initilized
            for i in range(len(row_idx)-1):
                num_nonzero_ele = row_idx[i+1] - row_idx[i]   # How many non-zero values in this row
                
                # To know which element of vector to product, we need to use "col" 
                # and "cursor".
                for j in range(num_nonzero_ele):    
                    result[i] += val[cursor+j]*vec[col[cursor+j]-1]
                cursor += num_nonzero_ele
            return result
        
        # If it's symmetric, (Ax + (A^T)x) - diag*x = the result of the complete matrix
        # because the symmetric matrix just contains the elements of the lower triangle
        if self.symmetric == True:                                                  
            CSR_result = CS_multi(self.val_CSR, self.col, self.row_idx, vec)            
            CSC_result = CS_multi(self.val_CSC, self.row, self.col_idx, vec)
            result = [a + b for a, b in zip(CSR_result, CSC_result)]    # Sum the 2 vectors
            idx = 0
            result1 = []
            for a,b in zip(result, self.diag):         # Subtract the diagonal product with vector
                result1.append(a - b * vec[idx])       
                idx += 1
            return result1
        else:        
            result = CS_multi(self.val_CSR, self.col, self.row_idx, vec)  
            return result    # if not symmetric, directly do matrix-vecotr product
        
'''
x_nonsym = [1] * 1030
x_sym = [1] * 5357
x_test = [3,4,5]
A = CSR('aaa.mtx')
b = A.mat_multi(x_test)
'''