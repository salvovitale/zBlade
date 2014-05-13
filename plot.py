""" @package zBlade
In this module it is implemented the class Plot which plot the parametric curves    
"""

from matplotlib.pylab import *
# from bspline import Bspline

# class Param_plot:
#     def __init__(self, color = '-r', legend = 'Insert a legend'):
#         self._c = color
#         self._l = legend
        


class Plot:
    
    def __init__(self, curve, ccp = 'bo', ccpoly = '-r', ccurve = '-b', lcp = 'Control Point', lcpoly = 'Control Polygon ', lcurve = 'Parametric Curve', xlabel = 'x', ylabel = 'y', discr = 100):
        self._curve = curve
        self._discr = discr
        self._x, self._y = self._generate()
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._ccp = ccp
        self._ccpoly = ccpoly
        self._ccurve = ccurve
        self._lcp = lcp
        self._lcpoly = lcpoly
        self._lcurve = lcurve
        if curve.__class__.__name__ == 'Bspline' :
             self._xU, self._yU = self._generate_knots()
            
            
        
    def __call__(self):
        plot(self._curve.get_x(), self._curve.get_y(), self._ccp)
        plot(self._curve.get_x(), self._curve.get_y(), self._ccpoly)
        plot(self._x, self._y, self._ccurve)
#         if self._curve.__class__.__name__ == 'Bspline' :
#             plot(self._xU, self._yU, 'ro') 
        legend([self._lcp, self._lcpoly, self._lcurve])
        xlabel(self._xlabel)
        ylabel(self._ylabel)
        
                    
    def _generate(self):
        t = np.linspace(0.0, 1.0, self._discr)
        x = np.ones(self._discr)
        y = np.ones(self._discr)
        for i in xrange(len(t)):
            x[i], y[i], nul1 = self._curve(t[i])
        return x, y 
       
    def _generate_knots(self):
        U = self._curve.get_int_U()
        n = len(U)
        x = np.ones(n)
        y = np.ones(n)
        for i in xrange(n):
            x[i], y[i], nul1 = self._curve(U[i])
        return x, y 