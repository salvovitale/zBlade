"""
Created on May 13, 2014

@author: salvovitale

In this module it is implemented the a Geometrical modeller for a convergent-divergent Nozzle    


"""

from matplotlib.pylab import *
from point import Point
from bspline import Bspline

class ConvDef(object):
    """
    we define the basic dimensions for the Converging part of the nozzle
    The convergin nozzle is dimensionless on the Throat Area
    
    """


    def __init__(self, l = 2.0, A = 2.0, m = 0.0, ncp = 6, type = 1):
        """
        Constructor
        
        l-> lenght
        A-> inlet Area
        m-> inlet derivative
        ncp-> number of Control Point
        type -> stretching function
        """
        
        self._l = l
        self._A = A
        self._m = m
        self._ncp = ncp
        self._type = type
        self._xp = self._spacing()
        self._Pig = self._init_cp()
         
        self._P = self._cp()
 
 
    def _spacing(self):
        ncp, l, type =  self._ncp, self._l, self._type
        if type == 1:
            xp = equi_spaced(-l, 0.0, ncp)
            return xp
        elif type == 2:
            xp = cos_spaced(-l, 0.0, ncp)
            return xp
        elif type == 3:
            xp = sin_spaced(-l, 0.0, ncp)
            return xp      
    

    def _init_cp(self):
        
        xp, m, A = self._xp, self._m, self._A
        dy1 = (xp[1] -xp[0])*m
        p01 = Point(xp[0], A)
        p11 = Point(xp[1], A + dy1)
        pn1 = Point(0.0, 1.0)
        pn_11 = Point(xp[-2], 1.0)
        return [p01, p11, pn_11, pn1]    
    

    def _interpol_init_guess(self):
        Pig = self._Pig
        n = len(Pig)
        x = np.ones(n)
        y = np.ones(n)
        for i in xrange(n):
            x[i], y[i], nul1, nul2 = Pig[i].split()
        from scipy import interpolate
        f = interpolate.lagrange(x, y)
        self._yp = f(self._xp)    
    
    
    def _cp(self):
        
        """ 
                    
        """
        self._interpol_init_guess()
        xp, yp = self._xp, self._yp    
        pct = []
        for i in xrange(len(xp)):
            pct.append(Point(xp[i], yp[i]))
 
        return pct
    
    def getCP(self):
        return self._P
    
    def getLenght(self): 
        return self._l
    
    def getArea(self): 
        return self._A        
    
    def getM(self): 
        return self._m
    
    def getNCP(self):
        return self._ncp

      
    

class DivDef(ConvDef):
    
    
    def _spacing(self):
        ncp, l, type =  self._ncp, self._l, self._type
        if type == 1:
            xp = equi_spaced(0.0, l, ncp)
            return xp
        elif type == 2:
            xp = cos_spaced(0.0, l, ncp)
            return xp
        elif type == 3:
            xp = sin_spaced(0.0, l, ncp)
            return xp        
    
    def _init_cp(self):
        
        xp, m, A = self._xp, self._m, self._A
        p01 = Point(xp[0], 1.0)
        p11 = Point(xp[1], 1.0)
        pn1 = Point(xp[-1], A)
        dy = (xp[-2] -xp[-1])*m
        pn_11 = Point(xp[-2], A + dy )
        return [p01, p11, pn_11, pn1]            

 


class ConvDiv:
    
    def __init__(self, cdDef, p = 3):
        
        self._cdDef = cdDef
        self._P = cdDef.getCP()
        self._p = p
        self._convDiv = Bspline(self._P, p = p)
        
    def __call__(self, u):
        return self._convDiv(u)
    
    def plot(self):
        self._convDiv.plot()
    
    def getCP(self):
        return self._P       







 
 
def equi_spaced(xin, xout, ndisc):
    
    dx = float(xout - xin)/(ndisc - 1)
    x = np.ones(ndisc)
    for i in xrange(ndisc):
                x[i] = xin +  i*dx 
    return x

def cos_spaced(xin, xout, ndisc):
    """
    this accumulate in the outlet
    """
    dalpha = (90.0)/(ndisc - 1)
    x = np.ones(ndisc)
    for i in xrange(ndisc):
                x[i] = xin + (xout - xin)*sin((i*dalpha)*pi/180.00) 
    return x

def sin_spaced(xin, xout, ndisc):
    """
    this accumulate in the inlet
    """
    dalpha = (90.0)/(ndisc - 1)
    x = np.ones(ndisc)
    for i in xrange(ndisc):
                x[ndisc-i-1] = xout - (xout - xin)*sin((i*dalpha)*pi/180.00) 
    return x        
               