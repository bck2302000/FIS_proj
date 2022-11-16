# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 09:11:04 2021

@author: Bo
"""
from Read_Matrix import CSR
import time
import matplotlib.pyplot as plt

class CG:
    def __init__(self, A, b, x0, threshold):
        self.A = A
        self.b = b
        self.x0 = x0
        self.r0 = [i - j for i, j in zip(b, A.mat_multi(x0))]
        self.p0 = self.r0
        self.threshold = threshold
        
    def run(self):
        alpha = []
        beta = []
        r = [self.r0]
        x = [self.x0]
        p = [self.p0]
        m = 0
        while 1:
            Ap = self.A.mat_multi(p[m])
            rm_dot = multiply(r[m],r[m], True)
            alpha.append(rm_dot / multiply(Ap, p[m], True))
            alpha_p = [alpha[m] * val for val in p[m]]
            #breakpoint()
            x.append([i + j for i, j in zip(x[m], alpha_p)])
            r.append([i - j for i, j in zip(r[m], [alpha[m] * val for val in Ap])])
            beta.append(multiply(r[m+1], r[m+1], True) / rm_dot)
            p.append([i + j for i, j in zip(r[m+1], [beta[m] * val for val in p[m]])])
            rho_data.append(L2_norm(r[m]))
            print('%s th iteration' %m)
            print(L2_norm(r[m]))
            #breakpoint()
            if L2_norm(r[m]) < self.threshold:
                m += 1
                break
            m += 1
            
        return m, L2_norm(r[m-1])
            
            
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

A = CSR('s3rmt3m3.mtx')  
#A = CSR('aaa.mtx')
x_star = [1] * A.shape[0]   #np.ones(A.shape[0])   
b = A.mat_multi(x_star)
x0 = [0] * A.shape[0]
threshold = 1e-8
rho_data = []
start = time.process_time()
CG_test = CG(A, b, x0, threshold)
num_iter, rho = CG_test.run()
end = time.process_time()
print('spent %s seconds' %(end - start))

x_axis = [i for i in range(1, num_iter + 1)]
plot(1, x_axis, rho_data, 'CG')