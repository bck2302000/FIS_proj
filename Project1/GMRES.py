# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 22:30:07 2021

Use GMRES to solve linear system Ax=b

@author: Bo
"""
from Read_Matrix import CSR

import time
import matplotlib.pyplot as plt
class GMRES:
    def __init__(self, matrixA: CSR, vectorb: list, m, threshold: float):
        self.A = matrixA
        self.b = vectorb
        self.threshold = threshold
        self.m = m
        self.precondition = 0                                          # 0 for without precondition, 1 for Jacobi, 2 for Gauss-Seidel
        if self.precondition == 2:
            self.GS_val, self.GS_col, self.GS_row_idx = self.A.Gauss_Seidel(self.A.FileName)
        
    def initial_guess(self, x0):
        self.x0 = x0
        self.r0 = [a - b for a,b in zip(self.b, self.A.mat_multi(x0))]          # self.b - self.A.mat_multi(x0)
        self.beta = L2_norm(self.r0)
# =============================================================================
#     def Gram_Schmidt(self, A, r0, m):
#         v = {}
#         for i in range(m):
#             v[i+1] = np.zeros_like(r0)
#         h = np.zeros((m+1, m))
#         v[1] = r0/np.linalg.norm(r0)
#         for j in range(1,m+1):
#             w = A.mat_multi(v[j])
#             for i in range(1,j+1):
#                 h[i][j] = sum(v[i] * w)
#                 w -= h[i][j]*v[i]
#             h[j+1][j] = np.linalg.norm(w)
#             v[j+1] = w/h[j+1][j]
#         return v, h
# 
# =============================================================================
    def run(self):
        v = [[0 for i in range(self.m + 1)] for j in range(len(self.r0))]        #np.zeros((len(self.r0), self.m+1))
        g = [0] * (self.m + 1) #np.zeros(self.m+1)
        g[0] = self.beta
        
        h = [[0 for i in range(self.m)] for j in range(self.m + 1)]            #np.zeros((self.m+1,self.m))
        #breakpoint()
        column_res(v, 0, [value / self.beta for value in self.r0])             #v[:,0] = self.r0/self.beta
        if self.precondition == 1 or self.precondition == 2:
            z = [[0 for i in range(self.m + 1)] for j in range(len(self.r0))]  #np.zeros((len(self.r0), self.m+1))
        c, s = {}, {}
        #self.A.Gauss_Seidel(A.FileName)
        global num_iter
        def GetKrylov(A, v):
            #breakpoint()
            if self.precondition == 1:
                column_res(z, j, multiply([1/val for val in A.diag], column(v, j)))        # z[:,j] = multiply(1/A.diag, v[:,j])                   
                w = A.mat_multi(column(z, j))                         # z[:,j])
                #breakpoint()
            elif self.precondition == 2:
                column_res(z, j, self.GS_forward_sub(self.GS_val, self.GS_col, self.GS_row_idx, column(v, j)))
                w = A.mat_multi(column(z, j))
                #breakpoint()
            else:    
                w = A.mat_multi(column(v, j))
                #breakpoint()
                
            for i in range(j+1):
                h[i][j] = multiply(w, column(v, i), scale = True)
                #breakpoint()
                #w -= h[i][j] * v[:, i]
                dummy = [h[i][j] * value for value in column(v, i)]
                w = [a - b for a,b in zip(w, dummy)]
                #w = [a - b for a,b in zip(w, [column(v,i)[idx]*h[i][j] for idx in range(len(column(v,i)))])]# h[i][j]*column(v, i)  #v[:, i]
            #breakpoint()
            h[j+1][j] = L2_norm(w)
            #breakpoint()
            column_res(v, j+1, [value / h[j+1][j] for value in w])                              #v[:, j+1] = w/h[j+1][j]
        
        for j in range(self.m):
            GetKrylov(self.A, v)
            if  full == True and j != 0: 
                if (self.precondition == 1 or self.precondition == 2):
                    ortho_data.append(abs(multiply(column(z, 0), column(z,j), scale = True)))                               #z[:,0], z[:, j], scale = True)))
                else:
                    ortho_data.append(abs(multiply(column(v, 0), column(v,j), scale = True)))
                
                
            for k in range(1,j+1):
                temp = h[k-1][j]
                h[k-1][j] = c[k-1] * h[k-1][j] + s[k-1] * h[k][j]
                h[k][j] = -s[k-1] * temp + c[k-1] * h[k][j]
                #breakpoint()
                
            c[j] = h[j][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
            s[j] = h[j+1][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
            h[j][j] = c[j] * h[j][j] + s[j] * h[j+1][j]
            h[j+1][j] = 0
            g[j+1] = -s[j] * g[j]
            g[j] = c[j] * g[j]
            #breakpoint()
            rho = abs(g[j+1]) / norm_r0
            #breakpoint()
            rho_data.append(rho)
            print('%s th iteration' %num_iter)
            Rm = [row[:j+1] for row in h[:j+1]]                              #for h[:j+1, :j+1]
            y_star = self.back_sub(Rm, g[:j+1])
            #breakpoint()
            if self.precondition == 1 or self.precondition == 2:
                vy = matmul([row[:j+1] for row in z], y_star)                             #z[:, :j+1]
               # breakpoint()
            else:
                vy = matmul([row[:j+1] for row in v], y_star)                             #v[:, :j+1]
            #breakpoint()
            xm = [a+b for a,b in zip(self.x0, vy)]
            error =  L2_norm([a - b for a, b in zip(x_star, xm)])
            error_data.append(error)
            if (rho < threshold): break
            num_iter += 1
            
        return xm, rho, error
    
    def back_sub(self, Rm, g):
        n = len(g)
        y = [0] * n  #np.zeros(n)
        for i in range(n-1, -1, -1):
            a = 0
            for j in range(i, n-1):
                a += Rm[i][j+1] * y[j+1]
            y[i] = (g[i] - a) / Rm[i][i]
        return y
    '''
    def GS_forward_sub(self, val, col, row_idx, v):
        n = len(v)
        z = np.zeros(n)
        for i in range(n):
            a = 0
            for j in range(i):
                a += M[i, j] * z[j]
            z[i] = (v[i] - a) / M[i, i]
        return z
    '''
    ### Use forward_sub to get the Z from MZ = v with the matrix is in CSR form
    def GS_forward_sub(self, val, col, row_idx, v):
        n = len(v)
        z = [0] * n  #np.zeros(n)
        for i in range(len(row_idx)-1):
            num_in_row = row_idx[i+1] - row_idx[i]
            a = 0
            for j in range(num_in_row-1):
                a += val[row_idx[i]+ j - 1] * z[col[row_idx[i] + j - 1] - 1]
            z[i] = (v[i] - a) / val[row_idx[i+1] - 2]
        return z
    
def plot(idx, x, y, title):
    plt.figure(idx)
    plt.plot(x, y)
    plt.yscale('log')
    plt.title(title)

def L2_norm(v):
    n = len(v)
    a = 0
    for i in range(n):
        a += v[i]**2
    return a**0.5

def multiply(v1, v2, scale = False):
    n = len(v1)
    a = [0] * n   #np.zeros(n)
    for i in range(n):
        a[i] = v1[i] * v2[i]
    #breakpoint()
    if scale == True:
        return sum(a)
    else:
        return a
    
def matmul(m1, v2):
    n = len(v2)
    num_row = len(m1)
    a = [0] * num_row    #np.zeros(num_row)
    #breakpoint()
    #if n == 1:
    #    return [m1[i] * v2 for i in range(num_row)]
    for i in range(num_row):
        for j in range(n):
            #breakpoint()
            a[i] += m1[i][j] * v2[j]
    return a

def column_res(matrix, idx, vec):
    for i, row in enumerate(matrix):
        row[idx] = vec[i]
# =============================================================================
#             breakpoint()
#             row[idx] = vec[idx]
# =============================================================================
            
def column(matrix, i):
    return [row[i] for row in matrix]        
    
A = CSR('orsirr_1.mtx')  
x_star = [1] * A.shape[0]   #np.ones(A.shape[0])   
b = A.mat_multi(x_star)
#b =  np.array([3.0, 2.0, 1.0])
x0 = [0] * A.shape[0]       #np.zeros(A.shape[0])    
#x0 = np.ones(3)
threshold = 1e-8

### Full GMRES
full = False


if full == True:
    num_iter = 1
    rho_data, ortho_data, error_data = [], [], []
    start = time.process_time()
    GMRES_test = GMRES(A, b, 600, threshold)
    GMRES_test.initial_guess(x0)
    norm_r0 = GMRES_test.beta
    xm, rho, error = GMRES_test.run()
    error = L2_norm([a - b for a, b in zip(xm, x_star)])
    end = time.process_time()
    x_axis = [i for i in range(1, num_iter + 1)]  #np.arange(1,num_iter+1)
    plot(1, x_axis, rho_data, 'relative residual with full GMRES' )
    plot(2, x_axis, error_data, 'actual residual with full GMRES' )
    plot(3, x_axis[1:], ortho_data, 'orthogonality with full GMRES')
    print('spent %s seconds' %(end - start))

### Restarted GMRES    
else:   
    m = 12
    num_iter = 1
    rho_data, error_data = [], []
    start = time.process_time()
    GMRES_test = GMRES(A, b, m, threshold)
    GMRES_test.initial_guess(x0)
    norm_r0 = GMRES_test.beta
    
    while 1:
        xm, rho, error = GMRES_test.run()
        if (rho < GMRES_test.threshold): break
        GMRES_test.initial_guess(xm)
        
    end = time.process_time()
    x_axis = [i for i in range(1, num_iter + 1)]  #np.arange(1,num_iter+1)
    plot(1, x_axis, rho_data, 'relative residual with restart number %s and Jacobi' %m)
    plot(2, x_axis, error_data, 'actual residual with restart number %s and Jacobi' %m)

    print('spent %s seconds' %(end - start))

