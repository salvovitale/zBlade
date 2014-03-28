""" @package zBlade
In this module it is implemented the a Geometrical modeller for a convergent-divergent Nozzle    
"""
# import sys
# utility_folder = '../utility/'
# sys.path.insert(0, utility_folder)
from matplotlib.pylab import *
import preproc_prof as pp
import scipy.optimize
from point import Point
from bezier import Bezier, curvi_abscissa
from bspline import Bspline



    
        
        
        
##################################################  NOOZLE GEO. MODELER ###############################################################        
    
         
class Conv(Bspline):
    def __init__(self, lc = 2.0, Ain = 2.0, min = -0.2, nc = 5, type = 1, p = 2):
        self._lc = lc
        self._Ain = Ain
        self._min = min
        self._nc = nc
        self._type = type
        self._Pc = self._pc()
        Bspline.__init__(self, self._Pc, p)
        
    def _pc(self):
        min, nc, lc, type = self._min, self._nc, self._lc, self._type
        """ 
            equispaced control point along x
        
        """
        if type == 1:
            dx = lc/(nc-1)
            dy1 = dx*min
            p01 = Point(-lc, self._Ain)
            p11 = p01 + Point(dx, dy1, 0.0, 0.0)
            pn1 = Point(0.0, 1.0)
            pn_11 = pn1 - Point(dx, 0.0, 0.0, 0.0)
            pint = [p01, p11, pn_11, pn1]
            n = len(pint)
            x = np.ones(n)
            y = np.ones(n)
            for i in xrange(n):
                x[i], y[i], nul1, nul2 = pint[i].split()
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
     
    
             




class Div(Bspline):
    
    def __init__(self, ld = 5.0, Aout = 3.0, mout = 0.1, nd = 5, type = 1, p = 2):
        self._ld = ld
        self._Aout = Aout
        self._mout = mout 
        self._nd = nd
        self._type = type
        Bspline.__init__(self, self._pd(), p)
      
    def _pd(self):
        mout, nd, ld, type = self._mout, self._nd,  self._ld, self._type
        if type == 1:
            dx = ld/(nd-1)
            p02 = Point(0, 1.0)
            p12 = p02 + Point(dx, 0.0, 0.0, 0.0)
            pn2 = Point( ld, self._Aout)
            dy2 = mout*dx
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
    
    def __init__(self, lc = 2.0, ld = 5.0, Ain = 2.0, Aout = 3.0, min = -0.2,  mout = 0.1, nc = 6, nd = 6, type = 1, p = 2):  
        """
        construct 2D - nozzle with 2 Bezier curve

        lc == length convergent part
        ld == length divergent part
        Ain == passage inlet area
        Aout == passage outlet areo
        
        all these quantities are adimensinal over Ath
        
        """      
        self._conv = Conv(lc, Ain, min, nc, type, p)
        self._div = Div(ld, Aout, mout, nd, type, p)

    def plot(self):
        self._conv.plot()
        self._div.plot()   
                

class BF_div(Div):
    
    def __init__(self, xy_p, ncp = 4, p = 3):
        self._xy_p = xy_p
        self._s = curvi_abscissa(xy_p)
        self._p = p
        self._init_curve(ncp)
        
        
    
    def _init_curve(self, ncp):
        _ld, _mout , _Aout = self._init_param()    
        Div.__init__(self, ld =_ld, Aout = _Aout, mout = _mout, nd = ncp, p = self._p)
        
    def _init_param(self):
        xy_p = self._xy_p
        l = len(xy_p[0,:])
        ld = xy_p[0,l-1] - xy_p[0,0]
        Aout = xy_p[1,l-1]
        mout = (xy_p[1,l-1] - xy_p[1,l-2])/(xy_p[0,l-1] - xy_p[0,l-2])
        return ld, mout, Aout    

    
    def __call__(self, A):
        P = self.get_P(A)
        self._fit = Bspline(P)
        return np.linalg.norm(self._err(), ord = 3.0)   
        
    
        
    def _err(self):
        xy_p, s =  self._xy_p, self._s
        x = np.ones(len(s))
        y = np.ones(len(s))
        err = np.ones(len(s))
        for i in xrange(len(s)):
            x[i], y[i], nul1  = Bspline.__call__(self, s[i])
            
            err[i] = abs(xy_p[0,i]- x[i]) + abs(xy_p[1,i]- y[i]) #it may be vectorized
            
#         err = abs(xy_p[0,:])
        return err
            
        
 
    def get_P(self, A):
        x, y, z, w = self.get_x(), self.get_y(), self.get_z(), self.get_w()
#         print x,y
        n = self.get_ncp()
        for i in xrange(n-2):
            y[i+1] = A[i]
        for i in xrange(n):
            w[i] = A[n-2+i]
                
        P = []
        for i in xrange(n):
            P.append(Point(x[i], y[i], z[i], w[i]))
        return P
     
    def get_A(self, P):
        n = self.get_ncp()
        x, y, z, w = np.ones(n), np.ones(n), np.ones(n), np.ones(n)
        for i in xrange(n):
            x[i], y[i], z[i], w[i] = P[i].split() 
        A = np.ones(2*n-2)
        for i in xrange(n-2):
            A[i] = y[i+1]
        for i in xrange(n):
            A[n-2+i] = w[i]    
        return A    
    
    def get_bounds(self):
        xy_p, n = self._xy_p, self.get_ncp()
        B  = np.ones((2*n-2 , 2), dtype = float)
        B[:,0] *= 0.1*xy_p[1,0]
        B[:,1] *= 1.01*xy_p[1,len(xy_p[0,:])-1]
        return B
    
    
    
        
           




if __name__=='__main__':            
    filename = 'div.dat'           
    xy_p = pp.read_xy_p(filename) 
    p = 7
    fit = BF_div(xy_p, ncp = 8 , p = p)
    Po = fit.get_cp()
#     init_guess = Bspline(Po)
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
    print P
    fit_curve = Bspline(P, p)
    print fit_curve.get_U()
    fit_curve.plot()
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
    
    err = [0.765107371068, 0.381679555444, 0.371869452821, 0.332929094068, 0.307328456426]
    ncp = [4, 5, 6, 7, 8]
    
#     axis([0, xy_p[0,len(xy_p[0,:])-1], 0,xy_p[0,len(xy_p[0,:])-1] ])
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
    
                
            