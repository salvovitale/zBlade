#from pylab import *
import sys
utility_folder = '../../utility/'
sys.path.insert(0, utility_folder)
from matplotlib.pylab import *
from math import factorial
import util as ut

import scipy.optimize
# import numpy as np


def distance_2p(xy1,xy2):
    return  np.sqrt((xy1[0] - xy2[0])**2 + (xy1[1] - xy2[1])**2)



def read_xy_p (filename):
    data = ut.read_data(filename)
    return np.array([[row[0] for row in data], [row[1] for row in data]])


_unzip = lambda zipped: zip(*zipped) # unzip a list of tuples

def _C(n, k):
    # binomial coefficient == n! / (i!(n - i)!)
    return factorial(n) / (factorial(k) * factorial(n - k))



def curvi_abscissa(xy_p):
    dist = np.ones(len(xy_p[0,:])-1)
    t = np.ones(len(xy_p[0,:]))
    t[0] = 0.0
    for i in xrange(len(dist)):
        dist[i]= distance_2p(xy_p[:,i],xy_p[:,i+1])
#             print xy_p[:,i],xy_p[:,i+1]
        t[i+1] = t[i] + dist[i]
    t /= t[len(t)-1]
    return t

# class Point
class Point:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0 , w = 1.0 ):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        
    def distance(self, other):
        return sqrt((other.x-self.x)*(other.x-self.x)+(other.y-self.y)*(other.y-self.y)+(other.z-self.z)*(other.z-self.z))

    def length(self):
        return self.distance(Point(0.0, 0.0, 0.0))

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z + other.z)

    def __mul__(self, c):
        return Point(c*self.x, c*self.y, c*self.z)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other):
        return not (self == other)
    
    def towards(self, target, t):
        return Point((1.0-t)*self.x+t*target.x, (1.0-t)*self.y+t*target.y, (1.0-t)*self.z+t*target.z)
    
    def halfway(self, target):
        return Point((self.x+target.x).div2(), (self.y+target.y).div2(), (self.z+target.z).div2())

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
        return self.x, self.y, self.z, self.w
    
    
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
    

class Plot:
    
    def __init__(self, curve, ccp = 'bo', ccpoly = '-b', ccurve = '-r', lccp = 'Control Point', lccpoly = 'Control Polygon ', lccurve = 'Parametric Curve', xlabel = 'x', ylabel = 'y', discr = 100):
        self._curve = curve
        self._discr = discr
        self._x, self._y = self._generate()
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._ccp = ccp
        self._ccpoly = ccpoly
        self._ccurve = ccurve
        self._lccp = lccp
        self._lccpoly = lccpoly
        self._lccurve = lccurve
        
    def __call__(self):
        plot(self._curve.get_x(), self._curve.get_y(), self._ccp)
        plot(self._curve.get_x(), self._curve.get_y(), self._ccpoly)
        plot(self._x, self._y, self._ccurve)
        legend([self._lccp, self._lccpoly, self._lccurve])
        xlabel(self._xlabel)
        ylabel(self._ylabel)
        
                    
    def _generate(self):
        t = np.linspace(0.0, 1.0, self._discr)
        x = np.ones(self._discr)
        y = np.ones(self._discr)
        for i in xrange(len(t)):
            x[i], y[i], nul1 = self._curve(t[i])
        return x, y    

    
class Bezier:
    def __init__(self, P, discr = 100):
        """
        construct bezier curve

        P == list of control points
        """
        self._n = len(P) # control point iterator
        self._P = P
        self._X, self._Y, self._Z, self._W = self.sep()  
        self._bc = self._bn()
        
        
    def __call__(self, t):
        X, Y, Z, W, n, bc = self._X, self._Y, self._Z, self._W, self._n, self._bc
        ber = self._bernstein(t)
        x = sum(X*bc*W*ber)
        y = sum(Y*bc*W*ber)
        z = sum(Z*bc*W*ber)
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
    
    
        
        
        
##################################################  NOOZLE GEO. MODELER ###############################################################        
    
         
class Conv(Bezier):
    def __init__(self, lc = 2.0, Ain = 2.0, m1 = -0.2, m2 = 0.0, nc = 4, type = 1):
        self._lc = lc
        self._Ain = Ain
        self._m1 = m1
        self._m2 = m2
        self._nc = nc
        self._type = type
        self._Pc = self._pc()
        Bezier.__init__(self, self._Pc)
        
    def _pc(self):
        m1, m2, nc, lc, type = self._m1, self._m2, self._nc, self._lc, self._type
        """ 
            equispaced control point along x
        
        """
        if type == 1:
            dx = lc/(nc-1)
            dy1 = dx*m1
            p01 = Point(-lc, self._Ain)
            p11 = p01 + Point(dx, dy1, 0.0, 0.0)
            pn1 = Point(0.0, 1.0)
            dy2 = dx*m2
            pn_11 = pn1 - Point(dx,dy2, 0.0, 0.0)
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
     
    
             




class Div(Bezier):
    
    def __init__(self, ld = 5.0, Aout = 3.0, m2 = 0.0, m3 = 0.1, nd = 4, type = 1):
        self._ld = ld
        self._Aout = Aout
        self._m2 = m2
        self._m3 = m3
        self._nd = nd
        self._type = type
        Bezier.__init__(self, self._pd())
      
    def _pd(self):
        m2, m3, nd, ld, type = self._m2, self._m3, self._nd,  self._ld, self._type
        if type == 1:
            dx = ld/(nd-1)
            p02 = Point(0, 1.0)
            dy1 = m2*dx
            p12 = p02 + Point(dx, dy1, 0.0, 0.0)
            pn2 = Point( ld, self._Aout)
            dy2 = m3*dx
            pn_12 = pn2 - Point(dx, dy2, 0.0, 0.0)
            pint =  [p02, p12, pn_12, pn2]
            n = len(pint)
            x = np.ones(n)
            y = np.ones(n)
            for i in xrange(n):
                x[i], y[i], nul1, nul2 = pint[i].split()
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
         
# 
# class Bezier_Ar(Bezier):
#     def _sep(self):
#         P = self._P
#         x = P[0,:]
#         y = P[1,:]
#         return x, y        

class BF_div(Div):
    
    def __init__(self, xy_p, ncp = 4):
        self._xy_p = xy_p
        self._s = curvi_abscissa(xy_p)
        self._init_curve(ncp)
        
    
    def _init_curve(self, ncp):
        _ld, _m2, _m3 , _Aout = self._init_param()    
        Div.__init__(self, ld =_ld, Aout = _Aout, m2 = _m2, m3 = _m3, nd = ncp)
        
    def _init_param(self):
        xy_p = self._xy_p
        l = len(xy_p[0,:])
        ld = xy_p[0,l-1] - xy_p[0,0]
        Aout = xy_p[1,l-1]
        m2 = 0.0
        m3 = (xy_p[1,l-1] - xy_p[1,l-2])/(xy_p[0,l-1] - xy_p[0,l-2])
        return ld, m2, m3, Aout    

    
    def __call__(self, A):
        P = self.get_P(A)
        self._fit = Bezier(P)
        return np.linalg.norm(self._err(), ord = 2.0)   
        
    
        
    def _err(self):
        xy_p, s =  self._xy_p, self._s
        x = np.ones(len(s))
        y = np.ones(len(s))
        err = np.ones(len(s))
        for i in xrange(len(s)):
            x[i], y[i], nul1  = Bezier.__call__(self, s[i])
            
            err[i] = abs(xy_p[0,i]- x[i]) + abs(xy_p[1,i]- y[i]) #it may be vectorized
            
#         err = abs(xy_p[0,:])
        return err
            
        
 
    def get_P(self, A):
        x, y, z, w = self.get_x(), self.get_y(), self.get_z(), self.get_w()
#         print x,y
        n = self.get_ncp()
        for i in xrange(n-4):
            y[i+2] = A[i]
#         for i in xrange(n):
#             w[i] = A[n-2+i]
                
        P = []
        for i in xrange(n):
            P.append(Point(x[i], y[i], z[i]))
        return P
     
    def get_A(self, P):
        n = self.get_ncp()
        x, y, z, w = np.ones(n), np.ones(n), np.ones(n), np.ones(n)
        for i in xrange(n):
            x[i], y[i], z[i], w[i] = P[i].split() 
        A = np.ones(n-3)
        for i in xrange(n-4):
            A[i] = y[i+2]
#         for i in xrange(n):
#             A[n-2+i] = w[i]    
        return A    
    
    def get_bounds(self):
        xy_p, n = self._xy_p, self.get_ncp()
        B  = np.ones((n-3 , 2), dtype = float)
        B[:,0] *= -1.1*xy_p[1,0]
        B[:,1] *= 3*xy_p[1,len(xy_p[0,:])-1]
        return B
    
    
    
        
           




if __name__=='__main__':            
    filename = 'div_geometry/my_data.txt'           
    xy_p = read_xy_p(filename) 
    fit = BF_div(xy_p, 7)
    Po = fit.get_cp()
    init_guess = Bezier(Po)
#     pig = Plot(init_guess, ccp = 'ko', ccpoly = '-k', ccurve = '-y', lccp = 'IG Control Point', lccpoly = 'IG Control Polygon ', lccurve = ' IG Parametric Curve')
#     pig() 
    
#     print Po
    Ao = fit.get_A(Po)
#     print Ao
    B =fit.get_bounds()
#     print Po, Ao
#     print B
      
#     fit.get_curve().bplot()
    
    A = scipy.optimize.fmin_slsqp(fit, Ao, bounds = B, iter = 1000)
    print A
    P = fit.get_P(A)
    fit_curve = Bezier(P)
    fit_curve.plot()
    print fit_curve.get_cp()
#     
    plot(xy_p[0,:], xy_p[1,:], '-g')

#     p = Plot(fit_curve)
#     p()
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
#     axis([-3, 6, 0, 10 ])
#     d.plot()
#     e.plot()
    
#     err = [0.765107371068, 0.381679555444, 0.371869452821, 0.332929094068, 0.307328456426]
#     ncp = [4, 5, 6, 7, 8]
    
    axis([0, xy_p[0,len(xy_p[0,:])-1], 0,xy_p[0,len(xy_p[0,:])-1] ])
#     title('nozzle modeler')
#     savefig('nozzle_modeler.eps')
#     title('curve fitting')
#     savefig('8pt.eps')
#     plot(ncp, err, '-bo')
#     plot(ncp, err, '-r')
#     xlabel('Number of equispaced Control Points')
#     ylabel('Least square error')
#     axis([3, 9, 0, 1 ])
#     title('fitting error')
#     savefig('fitting_error.eps')
#     axis('equal')
#     t = 0.5
    
#     curve = np.array(c(ti) for ti in t)
#     for ti in t:
#     print c.bezier(100)    
    
#     plot(curve)
    show()
    
                
            