#from pylab import *
import sys
utility_folder = '../utility/'
sys.path.insert(0, utility_folder)
from matplotlib.pylab import *
from math import factorial
import util as ut
import preproc_prof as pp
import scipy.optimize
# import numpy as np

_unzip = lambda zipped: zip(*zipped) # unzip a list of tuples

def _C(n, k):
    # binomial coefficient == n! / (i!(n - i)!)
    return factorial(n) / (factorial(k) * factorial(n - k))

# class Point
class Point:
    def __init__(self, x = 0.0, y = 0.0):
        self.x = x
        self.y = y
    
    def distance(self, other):
        return sqrt((other.x-self.x)*(other.x-self.x)+(other.y-self.y)*(other.y-self.y))

    def length(self):
        return self.distance(Point(0.0, 0.0))

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y)

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y)

    def __mul__(self, c):
        return Point(c*self.x, c*self.y)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return not (self == other)
    
    def towards(self, target, t):
        return Point((1.0-t)*self.x+t*target.x, (1.0-t)*self.y+t*target.y)
    
    def halfway(self, target):
        return Point((self.x+target.x).div2(), (self.y+target.y).div2())

    def compare_lex(self, other):
        if self.x < other.x:
            return -1
        if self.x > other.x:
            return 1
        if self.y < other.y:
            return -1
        if self.y > other.y:
            return 1
        return 0

    def less_lex(self, p):
        return self.compare_lex(p) < 0

    def less_eq_lex(self, p):
        return self.compare_lex(p) <= 0
    
    def __repr__(self):
        return "Point(%s, %s)" % (self.x, self.y)      
    
    def split(self):
        return self.x, self.y
    
    
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
    def __init__(self, P, discr = 100):
        """
        construct bezier curve

        P == list of control points
        """
        self._n = len(P) # control point iterator
        self._P = P
        self._X, self._Y = self.sep()  
        self._bc = self._bn()
        self._discr = discr
        self._x, self._y = self._generate()
        
    def sep(self):
        n, P = self._n, self._P
        x = np.ones(n)
        y = np.ones(n)
        for i in xrange(n):
            x[i], y[i] = P[i].split()
             
        return x, y    
        
    
    
    def __call__(self, t):
        X, Y, n, bc = self._X, self._Y, self._n, self._bc
        berns = np.ones(n)
        for i in xrange(n):
            berns[i] = t**i*(1-t)**(n-1-i)
        
        x = sum(X*bc*berns)
        y = sum(Y*bc*berns)
        return x, y     
        
    def _generate(self):
        t = np.linspace(0.0, 1.0, self._discr)
        x = np.ones(self._discr)
        y = np.ones(self._discr)
        for i in xrange(len(t)):
            x[i], y[i] = self(t[i])
        return x, y
    
    def _bn(self):
        n = self._n
        bc = np.ones(n)
        for i in  xrange(n):
            bc[i] = _C(n-1,i)
        return bc 
       
    def get_p(self):
        return self._P
    
    def get_ncp(self):
        return self._n
        
    def get_x(self):
        return self._X
    
    
    def get_y(self):
        return self._Y
    
    def bplot(self):    
#         plot(self._X, self._Y)
        plot(self._X, self._Y, 'bo')
        plot(self._x, self._y)
        
        
        
##################################################  NOOZLE GEO. MODELER ###############################################################        
    
         
class Conv:
    def __init__(self, lc = 2.0, Ain = 2.0, m1 = -0.2, m2 = 0.0, nc = 4, type = 1):
        self._lc = lc
        self._Ain = Ain
        self._m1 = m1
        self._m2 = m2
        self._nc = nc
        self._type = type
        self._Pc = self._pc()
        self._curve = Bezier(self._Pc)
        
    def _pc(self):
        m1, m2, nc, lc, type = self._m1, self._m2, self._nc, self._lc, self._type
        """ 
            equispaced control point along x
        
        """
        if type == 1:
            dx = lc/(nc-1)
            dy1 = dx*m1
            p01 = Point(-lc, self._Ain)
            p11 = p01 + Point(dx, dy1)
            pn1 = Point(0.0, 1.0)
            dy2 = dx*m2
            pn_11 = pn1 - Point(dx,dy2)
            pint = [p01, p11, pn_11, pn1]
            n = len(pint)
            x = np.ones(n)
            y = np.ones(n)
            for i in xrange(n):
                x[i], y[i] = pint[i].split()
            xp = np.ones(nc)
            for i in xrange(nc):
                xp[i] = -lc +  i*dx     
            
            
        from scipy import interpolate
        f = interpolate.lagrange(x, y)
        yp = f(xp)
#         print x, y
#         print xp, yp
            
#         f = interpolate.interp1d(x, y)
        pct = []
        for i in xrange(nc):
            pct.append(Point(xp[i], yp[i]))
#         pct.append(pn1)  
        return pct        
     
    def plot(self):
        self._curve.bplot()
             




class Div:
    
    def __init__(self, ld = 5.0, Aout = 3.0, m2 = 0.0, m3 = 0.1, nd = 4, type = 1):
        self._ld = ld
        self._Aout = Aout
        self._m2 = m2
        self._m3 = m3
        self._nd = nd
        self._type = type
        self._Pd = self._pd()
        self._curve  = Bezier(self._Pd)
        
    def __call__(self,t):
         return self._curve(t)   
        
    def _pd(self):
        m2, m3, nd, ld, type = self._m2, self._m3, self._nd,  self._ld, self._type
        if type == 1:
            dx = ld/(nd-1)
            p02 = Point(0, 1.0)
            dy1 = m2*dx
            p12 = p02 + Point(dx, dy1)
            pn2 = Point( ld, self._Aout)
            dy2 = m3*dx
            pn_12 = pn2 - Point(dx, dy2)
            pint =  [p02, p12, pn_12, pn2]
            n = len(pint)
            x = np.ones(n)
            y = np.ones(n)
            for i in xrange(n):
                x[i], y[i] = pint[i].split()
            xp = np.ones(nd)
            for i in xrange(nd):
                xp[i] =  i*dx     
        
        
        from scipy import interpolate
        f = interpolate.lagrange(x, y)
        yp = f(xp)
#         print x, y
#         print xp, yp
            
#         f = interpolate.interp1d(x, y)
        pct = []
        for i in xrange(nd):
            pct.append(Point(xp[i], yp[i]))
#         pct.append(pn1)  
        return pct
        
    def plot(self):
        self._curve.bplot()  
    
    def get_curve(self):
        return self._curve    
    
    def get_cp(self):
        return self._curve.get_p()
    
    def get_ncp(self):
        return self._curve.get_ncp() 
    
    def get_x(self):
        return self._curve.get_x()
    
    
    def get_y(self):
        return self._curve.get_y()




     
class Nozzle:
    
    def __init__(self, lc = 2.0, ld = 5.0, Ain = 2.0, Aout = 3.0, m1 = -0.2, m2 = 0.0, m3 = 0.1, nc = 4, nd = 4, type = 1):  
        """
        construct 2D - nozzle with 2 Bezier curve

        lc == length convergent part
        ld == length divergent part
        Ain == passage inlet area
        Aout == passage outlet areo
        
        all these quantities are adimensinal over Ath
        
        """      
        self._conv = Conv(lc, Ain, m1, m2, nc, type)
        self._div = Div(ld, Aout, m2, m3, nd, type)

    def plot(self):
        self._conv.plot()
        self._div.plot()   
         

class Bezier_Ar(Bezier):
    def _sep(self):
        P = self._P
        x = P[0,:]
        y = P[1,:]
        return x, y        

class BF_div:
    
    def __init__(self, xy_p, ncp = 4):
        self._xy_p = xy_p
        self._t = self.curvi_abscissa(xy_p)
        self._ncp = ncp
        self._curve = self._init_curve()
    
        
    def __call__(self, A):
        P = self.get_P(A)
        self._fit = Bezier(P)
        return np.linalg.norm(self._err(), ord = 2.0)
        
    def _err(self):
        curve, xy_p, t = self._curve, self._xy_p, self._t
        x = np.ones(len(t))
        y = np.ones(len(t))
        err = np.ones(len(t))
        for i in xrange(len(t)):
            x[i], y[i] = curve(t[i])
            
            err[i] = abs(xy_p[0,i]- x[i]) + abs(xy_p[1,i]- y[i]) #it may be vectorized
            
#         err = abs(xy_p[0,:])
        return err
            
        
    def _init_curve(self):
        _ld, _m2, _m3 , _Aout = self._init_param()    
        return  Div(ld =_ld, Aout = _Aout, m2 = _m2, m3 = _m3, nd = self._ncp)
        
    def _init_param(self):
        xy_p = self._xy_p
        l = len(xy_p[0,:])
        ld = xy_p[0,l-1] - xy_p[0,0]
        Aout = xy_p[1,l-1]
        m2 = 0.0
        m3 = (xy_p[1,l-1] - xy_p[1,l-2])/(xy_p[0,l-1] - xy_p[0,l-2])
        return ld, m2, m3, Aout



    def get_P(self, A):
        x, y = self._curve.get_x(), self._curve.get_y()
        n = self._curve.get_ncp()
        for i in xrange(n-2):
            y[i+1] = A[i]
        P = []
        for i in xrange(n):
            P.append(Point(x[i], y[i]))
        return P
    
    def get_A(self, P):
        n = self._curve.get_ncp()
        x, y = np.ones(len(P)), np.ones(len(P))
        for i in xrange(n):
            x[i], y[i] = P[i].split() 
        A = np.ones(n-2)
        for i in xrange(n-2):
            A[i] = y[i+1]
        return A    
    
    def get_bounds(self):
        xy_p, n = self._xy_p, self._curve.get_ncp()
        B  = np.ones((n-2 , 2), dtype = float)
        B[:,0] *= xy_p[1,0]
        B[:,1] *= xy_p[1,len(xy_p[0,:])-1]
        return B
    
    def get_curve(self):
        return self._curve
    
    
        
           
    @staticmethod
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
#     @staticmethod
#     def init_cfc( cfc_ig1 = 0.4, cfc_ig2 = 0.4):
#         N = np.ones(2)
#         N[0] *= cfc_ig1
#         N[1] *= cfc_ig2
#         return N
 
#     @staticmethod
#     def init_A(poly_gr = 5, A_ig = 0.8):
#         A = np.ones(poly_gr+1)
#         A *= A_ig
#         return A
# 
#     @staticmethod    
#     def bound_cfc( cfc_lb1 = 0.1, cfc_ub1 = 1.5, cfc_lb2 = 0.1 , cfc_ub2 = 1.4):
#         N = np.ones((2,2), dtype = float)
#         N[0,0] *= cfc_lb1
#         N[1,0] *= cfc_lb2
#         N[0,1] *= cfc_ub1
#         N[1,1] *= cfc_ub2
#         return N


if __name__=='__main__':            
    filename = 'div.dat'           
    xy_p = pp.read_xy_p(filename) 
    fit = BF_div(xy_p, 8)
    Po = fit.get_curve().get_cp()
    Ao = fit.get_A(Po)
    B =fit.get_bounds()
#     print Po, Ao
    print B
    plot(xy_p[0,:], xy_p[1,:] )  
    fit.get_curve().plot()
    
    A = scipy.optimize.fmin_slsqp(fit, Ao, bounds = B, iter = 1000)
    print A
    P = fit.get_P(A)
    fit_curve = Bezier(P)
    fit_curve.bplot()
#     print BF_div.curvi_abscissa(xy_p) 
    """
        creare una Bezier che usa array e non punti ed il gioco e fatto
        usare la classe div per trovare i punti di controllo iniziali
    """       
#     p0 = Point(0.0, 0.0)
#     p1 = Point(0.5, 0.5)
#     p2 = Point(1.0, 0.0)
#     p3 = Point(1.0, 1.0)
#     P= [p0, p1, p2, p3]
#     
#     c = Nozzle()
#     d = Nozzle(nc = 5, nd = 5)
#     e = Nozzle(nc = 6, nd = 6)
#     c.plot()
#     d.plot()
#     e.plot()
    axis('equal')
#     t = 0.5
    
#     curve = np.array(c(ti) for ti in t)
#     for ti in t:
#     print c.bezier(100)    
    
#     plot(curve)
    show()
    
                
            