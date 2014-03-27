""" @package zBlade
In this module it is implemented the class Bezier   
"""

import preproc_prof as pp
import numpy as np
from math import factorial
from plot import Plot


_unzip = lambda zipped: zip(*zipped) # unzip a list of tuples

def _C(n, k):
    # binomial coefficient == n! / (i!(n - i)!)
    return factorial(n) / (factorial(k) * factorial(n - k))



def curvi_abscissa(xy_p):
    dist = np.ones(len(xy_p[0,:])-1)
    t = np.ones(len(xy_p[0,:]))
    t[0] = 0.0
    for i in xrange(len(dist)):
        dist[i]= pp.distance_2p(xy_p[:,i],xy_p[:,i+1])
#             print xy_p[:,i],xy_p[:,i+1]
        t[i+1] = t[i] + dist[i]
    t /= t[len(t)-1]
    return t

# class Point

    
class Lagrange:

    def __init__(self, P, t):
        """
        construct Lagrange curve function

        P == list of control points
        t == list of time points
        len(P) == len(t)
        """
        n = len(t)
        assert len(P) == n # same number of time and control points
        self.X, self.Y = _unzip(P) # points in X, Y components
        self._n = xrange(n) # time/control point iterator
        self.t = t

    def __call__(self, t_):
        """
        return point on Langrange function at t_
        """
        X, Y, t, _n = self.X, self.Y, self.t, self._n
        x, y = 0, 0 # initial x and y return values
        for i in _n:
            p_i = 1 # initial lagrange polynomial value
            for j in _n:
                # if i != j: update lagrange polynomial
                if i != j: p_i *= (t_ - t[j]) / (t[i] - t[j])
            # mult ith control point by ith lagrange polynomial
            # (ith control point maps to ith time point)
            x += X[i] * p_i
            y += Y[i] * p_i
        return x, y    
    
  
class Bezier:
    def __init__(self, P):
        """
        construct  rational bezier curve

        P == list of control points
        """
        self._n = len(P) # control point iterator
        self._P = P
        self._X, self._Y, self._Z, self._W = self.sep()  
        self._bc = self._bn()
        
        
    def __call__(self, t):
        X, Y, Z, W, n, bc = self._X, self._Y, self._Z, self._W, self._n, self._bc
        ber = self._bernstein(t)
        x = sum(X*bc*ber)
        y = sum(Y*bc*ber)
        z = sum(Z*bc*ber)
        w = sum(W*bc*ber)
        return x/w, y/w, z/w    
        
    def sep(self):
        n, P = self._n, self._P
        x = np.ones(n)
        y = np.ones(n)
        z = np.ones(n)
        w = np.ones(n)
        for i in xrange(n):
            x[i], y[i], z[i], w[i] = P[i].split()
             
        return x, y, z, w           
    
    
    
    def _bn(self):
        n = self._n
        bc = np.ones(n)
        for i in  xrange(n):
            bc[i] = _C(n-1,i)
        return bc 
       
    def _bernstein(self, t):
        n = self._n
        ber = np.ones(n)
        for i in xrange(n):
            ber[i] = t**i*(1-t)**(n-1-i)
        return ber    
           
    def get_cp(self):
        return self._P
    
    def get_ncp(self):
        return self._n
        
    def get_x(self):
        return self._X
    
    
    def get_y(self):
        return self._Y
    
    def get_z(self):
        return self._Z
    
    
    def get_w(self):
        return self._W
    
    def plot(self):
        a = Plot(self)
        a()
    