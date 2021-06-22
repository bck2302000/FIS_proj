# -*- coding: utf-8 -*-
"""
Created on Sun Jun 20 17:38:46 2021

@author: Bo
"""
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 22:30:07 2021

Use GMRES to solve linear system Ax=b

@author: Bo
"""
from Read_Matrix_np import CSR
import numpy as np
import time
import matplotlib.pyplot as plt
class GMRES:
    def __init__(self, matrixA: CSR, vectorb: np.array, m, threshold: float):
        self.A = matrixA
        self.b = vectorb
        self.threshold = threshold
        self.m = m
        self.precondition = 2                                          # 0 for without precondition, 1 for Jacobi, 2 for Gauss-Seidel
        if self.precondition == 2:
            self.GS_val, self.GS_col, self.GS_row_idx = self.A.Gauss_Seidel(self.A.FileName)
        
    def initial_guess(self, x0):
        self.x0 = x0
        self.r0 = self.b - self.A.mat_multi(x0)
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
        v = np.zeros((len(self.r0), self.m+1))
        g = np.zeros(self.m+1)
        g[0] = self.beta
        h = np.zeros((self.m+1,self.m))
        v[:, 0] = self.r0/self.beta
        if self.precondition == 1 or self.precondition == 2:
            z = np.zeros((len(self.r0), self.m+1))
        c, s = {}, {}
        #self.A.Gauss_Seidel(A.FileName)
        global num_iter
        def GetKrylov(A, v):
            if self.precondition == 1:
                z[:,j] = multiply(1/A.diag, v[:, j])
                w = A.mat_multi(z[:,j])
            elif self.precondition == 2:
                z[:,j] = self.GS_forward_sub(self.GS_val, self.GS_col, self.GS_row_idx, v[:, j])
                w = A.mat_multi(z[:,j])
            else:    
                w = A.mat_multi(v[:,j])
                
            for i in range(j+1):
                h[i][j] = multiply(w, v[:, i], scale = True)
                #breakpoint()
                w -= h[i][j]*v[:, i]
            h[j+1][j] = L2_norm(w)
            breakpoint()
            v[:, j+1] = w/h[j+1][j]
        
        for j in range(self.m):
            GetKrylov(self.A, v)
            
            if  full == True and j != 0: 
                if (self.precondition == 1 or self.precondition == 2):
                    ortho_data.append(abs(multiply(z[:,0], z[:, j], scale = True)))
                else:
                    ortho_data.append(abs(multiply(v[:,0], v[:, j], scale = True)))
                
                
            for k in range(1,j+1):
                temp = h[k-1][j]
                h[k-1][j] = c[k-1] * h[k-1][j] + s[k-1] * h[k][j]
                h[k][j] = -s[k-1] * temp + c[k-1] * h[k][j]
                breakpoint()
                
            c[j] = h[j][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
            s[j] = h[j+1][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
            h[j][j] = c[j] * h[j][j] + s[j] * h[j+1][j]
            h[j+1][j] = 0
            g[j+1] = -s[j] * g[j]
            g[j] = c[j] * g[j]
            breakpoint()
            rho = abs(g[j+1]) / norm_r0
            rho_data.append(rho)
            print('%s th iteration' %num_iter)
            Rm = h[:j+1, :j+1]
            y_star = self.back_sub(Rm, g[:j+1])
            if self.precondition == 1 or self.precondition == 2:
                vy = matmul(z[:, :j+1], y_star)
               # breakpoint()
            else:
                vy = matmul(v[:, :j+1], y_star)
            xm = self.x0 + vy
            error =  L2_norm(x_star - xm)
            error_data.append(error)
            if (rho < threshold): break
            num_iter += 1
            
        return xm, rho, error
    
    def back_sub(self, Rm, g):
        n = len(g)
        y = np.zeros(n)
        for i in range(n-1, -1, -1):
            a = 0
            for j in range(i, n-1):
                a += Rm[i, j+1] * y[j+1]
            y[i] = (g[i] - a) / Rm[i, i]
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
        z = np.zeros(n)
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
    a = np.zeros(n)
    for i in range(n):
        a[i] = v1[i] * v2[i]
    if scale == True:
        return sum(a)
    else:
        return a
    
def matmul(m1, v2):
    n = len(v2)
    num_row = m1.shape[0]
    a = np.zeros(num_row)
    for i in range(num_row):
        for j in range(n):
            a[i] += m1[i , j] * v2[j]
    return a

A = CSR('orsirr_1.mtx')  
x_star = np.ones(A.shape[0])   
b = A.mat_multi(x_star)
#b =  np.array([3.0, 2.0, 1.0])
x0 = np.zeros(A.shape[0])    
#x0 = np.ones(3)
threshold = 1e-8

### Full GMRES
full = True


if full == True:
    num_iter = 1
    rho_data, ortho_data, error_data = [], [], []
    start = time.process_time()
    GMRES_test = GMRES(A, b, 600, threshold)
    GMRES_test.initial_guess(x0)
    norm_r0 = GMRES_test.beta
    xm, rho, error = GMRES_test.run()
    error = L2_norm(xm - x_star)
    end = time.process_time()
    x_axis = np.arange(1,num_iter+1)
    plot(1, x_axis, rho_data, 'relative residual with full GMRES' )
    plot(2, x_axis, error_data, 'actual residual with full GMRES' )
    plot(3, x_axis[1:], ortho_data, 'orthogonality with full GMRES')
    print('spent %s seconds' %(end - start))

### Restarted GMRES    
else:   
    m = 50
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
    x_axis = np.arange(1,num_iter+1)
    plot(1, x_axis, rho_data, 'relative residual with restart number %s and Gauss-Seidel' %m)
    plot(2, x_axis, error_data, 'actual residual with restart number %s and Gauss-Seidel' %m)

    print('spent %s seconds' %(end - start))

