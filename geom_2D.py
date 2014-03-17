#from pylab import *
from matplotlib.pylab import *
from math import factorial
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
        self._X, self._Y = self._sep()  
        self._bc = self._bn()
        self._discr = discr
        self._x, self._y = self._generate()
        
    def _sep(self):
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
        
    def get_x(self):
        return self._X
    
    
    def get_y(self):
        return self._Y
    
    def bplot(self):    
        plot(self._X, self._Y)
        plot(self._X, self._Y, 'bo')
        plot(self._x, self._y)
    
        
 
 
class Nozzle:
    
    def __init__(self, lc = 2.0, ld = 5.0, Ain = 2.0, Aout = 3.0, m1 = -0.2, m2 = 0.0, m3 = 0.1, nc = 3, nd = 3, type = 1):  
        """
        construct 2D - nozzle with 2 Bezier curve

        lc == length convergent part
        ld == length divergent part
        Ain == passage inlet area
        Aout == passage outlet areo
        
        all these quantities are adimensinal over Ath
        
        """      

        self._lc = lc
        self._ld = ld
        self._Ain = Ain
        self._Aout = Aout
        self._m1 = m1
        self._m2 = m2
        self._m3 = m3
        self._nc = nc
        self._nd = nd
        self._type
        self._Pc = self._pc()
        self._Pd = self._pd()
        self._conv = Bezier(self._Pc)
        self._div  = Bezier(self._Pd)
    
    
    def _pc(self):
        m1, m2, nc, lc, type = self._m1, self._m2, self._nc, self._lc, self.type
        if type == 1:
            dx = lc/nc
            dy1 = dx*m1
            p01 = Point(0.0, self._Ain)
            p11 = p01 + Point(dx, dy1)
            pn1 = Point(lc, 1.0)
            dy2 = dx*m2
            pn_11 = pn1 - Point(dx,dy2)
            pint = [p01, p11, pn_11, pn1]
            n = len(pint)
            x = np.ones(n)
            y = np.ones(n)
            for i in xrange(n):
                x[i], y[i] = pint[i].split()
            xp = np.ones[nc+1]
            for i in xrange(nc+1):
                xp[i] =     
        from scipy import interpolate
        f = interpolate.interp1d(x, y)
            
            
        f = interpolate.interp1d(x, y)  
        return [p01, p11, pn_11, pn1]
    
    
    def _pd(self):
        m2, m3, nd, lc, ld = self._m2, self._m3, self._nd, self._lc, self._ld
        dx = ld/nd
        p02 = Point(lc, 1.0)
        dy1 = m2*dx
        p12 = p02 + Point(dx, dy1)
        pn2 = Point(lc + ld, self._Aout)
        dy2 = m3*dx
        pn_12 = pn2 - Point(dx, dy2)
        
        return [p02, p12, pn_12, pn2]
        
    def plot(self):
        self._conv.bplot()
        self._div.bplot()    
if __name__=='__main__':            
            
            
#     p0 = Point(0.0, 0.0)
#     p1 = Point(0.5, 0.5)
#     p2 = Point(1.0, 0.0)
#     p3 = Point(1.0, 1.0)
#     P= [p0, p1, p2, p3]
#     
    c = Nozzle()
    
    c.plot()
    axis('equal')
#     t = 0.5
    
#     curve = np.array(c(ti) for ti in t)
#     for ti in t:
#     print c.bezier(100)    
    
#     plot(curve)
    show()
    
                
            