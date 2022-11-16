# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 23:24:12 2021

@author: Bo
"""

# Do an operation on the vector, dependent on the users selection
def vec_op(A: list, val: int, ope: str):
    res = [0] * len(A)
    for i in range(len(A)):
        res[i] = A[i]
    if ope == '+':
        for i in range(len(A)):
            res[i] += val
    elif ope == '-':
        for i in range(len(A)):
            res[i] -= val
    elif ope == '*':
        for i in range(len(A)):
            res[i] *= val
    elif ope == '/':
        for i in range(len(A)):
            res[i] /= val
    else:
        print('Wrong operator!')
    return res

# Access the L2_norm of a vector
def L2_norm(v):
    n = len(v)
    a = 0
    for i in range(n):
        a += v[i]**2
    return a**0.5

# multiply 2 vectors. Whether the result is scalar or a vector depends on the user.
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
  
# substract 2 vectors
def vec_subtract(v1, v2):
    n = len(v1)
    a = [0] * n
    for i in range(n):
        a[i] = v1[i] - v2[i]
    return a