# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 21:41:35 2021

@author: Bo
"""
#from Read_Matrix import CSR
import numpy as np


'''
def run(A, b, x0, m, threshold):
    r0 = b-np.matmul(A, x0)
    beta = np.linalg.norm(r0)
    v = np.zeros((len(r0), m+1))
    g = np.zeros(m+1)
    g[0] = beta
    h = np.zeros((m+1,m))
    v[:, 0] = r0/beta
    c, s = {}, {}
    def GetKrylov(A, v):
        w = np.matmul(A,v[:,j])
        for i in range(j+1):
            h[i][j] = np.matmul(w.transpose(), v[:, i])
            w -= h[i][j]*v[:, i]
        h[j+1][j] = np.linalg.norm(w)
        v[:, j+1] = w/h[j+1][j]
    
    for j in range(m):
        GetKrylov(A, v)
        for k in range(1,j+1):
            temp = h[k-1][j]
            h[k-1][j] = c[k-1] * h[k-1][j] + s[k-1] * h[k][j]
            h[k][j] = -s[k-1] * temp + c[k-1] * h[k][j]
        breakpoint()
        c[j] = h[j][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
        s[j] = h[j+1][j] / (h[j][j]**2 + h[j+1][j]**2)**0.5
        h[j][j] = c[j] * h[j][j] + s[j] * h[j+1][j]
        g[j+1] = -s[j] * g[j]
        g[j] = c[j] * g[j]
        
    Rm = h[:m, :]
    Rm_inv = np.linalg.inv(Rm)
    y_star = np.matmul(Rm_inv, g[:m])
    vy = np.zeros(len(r0))
    
    for i in range(len(self.r0)):
        for j in range(self.m):
            #vy[i] += v[j][i] * y_star[j]
            breakpoint()
    
    #vy = np.matmul(v[:, :self.m], y_star)
    for i in range(len(r0)):
        for j in range(m):
            vy[i] += v[i, j] * y_star[j]
    xm = x0 + vy
    return xm, np.linalg.norm(b - np.matmul(A,xm))
   
A = np.array([[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]])
x = np.ones(4)
b = np.matmul(A,x)
x0 = np.zeros(4)
m = 2
threshold = 0.01
xm, rho = run(A, b, x0, m, threshold)
while (rho > threshold):
    xm, rho = run(A, b, xm, 3, 0.01)
    print(rho)
'''       

#import numpy as np
import math

class GMRES_API(object):
    def __init__( self,
                  A_coefficient_matrix: np.array([], dtype = float ),
                  b_boundary_condition_vector: np.array([], dtype = float ),
                  maximum_number_of_basis_used: int,
                  threshold = 1.0e-16 ):

        self.A = A_coefficient_matrix
        self.b = b_boundary_condition_vector
        self.maximum_number_of_basis_used = maximum_number_of_basis_used
        self.threshold = threshold

    def initial_guess_input( self, x_input_vector_initial_guess: np.array([], dtype = float ) ):

        self.x = x_input_vector_initial_guess

        try:
            assert len( self.x ) == len( self.b )

        except Exception:

            print(" The input guess vector's size must equal to the system's size !\n")
            print(" The matrix system's size == ", len( self.b ))
            print(" Your input vector's size == ", len( self.x ))
            self.x = np.zeros( len( self.b ) ) 
            print(" Use default input guess vector = ", self.x, " instead of the incorrect vector you given !\n")


    def run( self ):

        n = len( self.A )
        m = self.maximum_number_of_basis_used

        r = self.b - np.dot(self.A , self.x)
        r_norm = np.linalg.norm( r )

        b_norm = np.linalg.norm( self.b )

        self.error = np.linalg.norm( r ) / b_norm
        self.e = [self.error]
        
        # initialize the 1D vectors 
        sn = np.zeros( m )
        cs = np.zeros( m )
        e1 = np.zeros( m + 1 )
        e1[0] = 1.0

        beta = r_norm * e1 
        # beta is the beta vector instead of the beta scalar

        H = np.zeros(( m+1, m+1 ))
        Q = np.zeros((   n, m+1 ))
        Q[:,0] = r / r_norm
        for k in range(m):

            ( H[0:k+2, k], Q[:, k+1] )    = __class__.arnoldi( self.A, Q, k)
            breakpoint()
            ( H[0:k+2, k], cs[k], sn[k] ) = __class__.apply_givens_rotation( H[0:k+2, k], cs, sn, k)
            
            # update the residual vector
            beta[ k+1 ] = -sn[k] * beta[k]
            beta[ k   ] =  cs[k] * beta[k]

            # calculate and save the errors
            self.error = abs(beta[k+1]) / b_norm
            self.e = np.append(self.e, self.error)

            if( self.error <= self.threshold):
                break


        # calculate the result
        #y = np.matmul( np.linalg.inv( H[0:k+1, 0:k+1]), beta[0:k+1] )
        #TODO Due to H[0:k+1, 0:k+1] being a upper tri-matrix, we can exploit this fact. 
        y = __class__.__back_substitution( H[0:k+1, 0:k+1], beta[0:k+1] )


        self.x = self.x + np.matmul( Q[:,0:k+1], y )
        breakpoint()
        self.final_residual_norm = np.linalg.norm( self.b - np.matmul( self.A, self.x ) )

        return self.x


    '''''''''''''''''''''''''''''''''''
    '        Arnoldi Function         '
    '''''''''''''''''''''''''''''''''''
    @staticmethod
    def arnoldi( A, Q, k ):
        h = np.zeros( k+2 )
        q = np.dot( A, Q[:,k] )
        for i in range ( k+1 ):
            h[i] = np.dot( q, Q[:,i])
            q = q - h[i] * Q[:, i]
        h[ k+1 ] = np.linalg.norm(q)
        q = q / h[ k+1 ]
        return h, q 

    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    '           Applying Givens Rotation to H col           '
    '''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    @staticmethod
    def apply_givens_rotation( h, cs, sn, k ):
        for i in range(1, k+1 ):
            temp   =  cs[i-1] * h[i-1] + sn[i-1] * h[i]
            h[i] = -sn[i-1] * h[i-1] + cs[i-1] * h[i]
            h[i-1]   = temp
            breakpoint()
        # update the next sin cos values for rotation
        cs_k, sn_k, h[k] = __class__.givens_rotation( h[k-1], h[k] )
        breakpoint()
        # eliminate H[ k+1, i ]
        h[k + 1] = 0.0

        return h, cs_k, sn_k

    ##----Calculate the Given rotation matrix----##
    # From "http://www.netlib.org/lapack/lawnspdf/lawn150.pdf"
    # The algorithm used by "Edward Anderson"
    @staticmethod
    def givens_rotation( v1, v2 ):
        if( v2 == 0.0 ):
            cs = np.sign(v1)
            sn = 0.0
            r = abs(v1)
        elif( v1 == 0.0 ):
            cs = 0.0
            sn = np.sign(v2)
            r = abs(v2)
        elif( abs(v1) > abs(v2) ):
            t = v2 / v1 
            u = np.sign(v1) * math.hypot( 1.0, t )  
            cs = 1.0 / u
            sn = t * cs
            r = v1 * u
        else:
            t = v1 / v2 
            u = np.sign(v2) * math.hypot( 1.0, t )  
            sn = 1.0 / u
            cs = t * sn
            r = v2 * u
        return cs, sn, r

    # From https://stackoverflow.com/questions/47551069/back-substitution-in-python
    @staticmethod
    def __back_substitution( A: np.ndarray, b: np.ndarray) -> np.ndarray:
        n = b.size
        if A[n-1, n-1] == 0.0:
            raise ValueError

        x = np.zeros_like(b)
        x[n-1] = b[n-1] / A[n-1, n-1]
        for i in range( n-2, -1, -1 ):
            bb = 0
            for j in range ( i+1, n ):
                bb += A[i, j] * x[j]
            x[i] = (b[i] - bb) / A[i, i]
        return x


    def final_residual_info_show( self ):
        print( "x  =", self.x, "residual_norm =  ", self.final_residual_norm ) 
    
def main():

    A_mat = np.array( [[1.00, 1.00, 1.00],
                       [1.00, 2.00, 1.00],
                       [0.00, 0.00, 3.00]] )

    b_mat = np.array( [3.0, 2.0, 1.0] )

    GMRES_test_itr2 = GMRES_API( A_mat, b_mat, 5, 0.01)

    x_mat = np.array( [1.0, 1.0, 1.0] )
    print("x  =", x_mat)

    # GMRES with restart, 2 iterations in each restart ( GMRES(2) )
    max_restart_counts = 100
    for restart_counter in range(max_restart_counts):
        GMRES_test_itr2.initial_guess_input( x_mat )

        x_mat = GMRES_test_itr2.run()
        print(restart_counter+1," : x  =", x_mat)

    xx = np.matmul( np.linalg.inv(A_mat), b_mat )
    print("ANS : xx =", xx) 


if __name__ == '__main__':
    main()  
'''
A = CSR('orsirr_1.mtx')     
b = A.mat_multi(np.ones(1030))
x0 = np.zeros(1030)    
GMRES_test = GMRES(A, b, 10, 1)
GMRES_test.initial_guess(x0)
rho = GMRES_test.beta
while (rho > GMRES_test.threshold):
    xm, rho = GMRES_test.run()
    print(rho)
    GMRES_test.initial_guess(xm)
'''

