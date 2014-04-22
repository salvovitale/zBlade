"""
Created on Apr 22, 2014

@author: salvovitale
"""




import sys
zblade_folder= '../'
sys.path.insert(0, zblade_folder)
from matplotlib.pylab import *
# from point import Point, Norm, Angle
# from line  import Line
from point import Point, Angle, Norm
from bspline import *
from line import Line
from bspline import Bspline
from math import sqrt, pi, cos, sin, tan

class Udist(object):
    """
    classdocs
    """


    def __init__(self, nintP, opt = 1, a = 0.0 , b = 1.0):
        
        """
        Constructor
        """
        self._opt = opt
        self._udistribution()
        self._nintP = nintP
        self._a = a
        self._b = b
    
    def __call__(self):
        return self._u
    
    
    def _udistribution(self):
        opt = self._opt
        if opt == 1:
            self._equispaced()
            
            
    def _equispaced(self):
        nintP, a, b = self._nintP, self._a, self._b
        d = (b-a)/(nintP+1)
        u = np.ones(nintP+2)
        for i in xrange(nintP+2):
            u.append(a + d*i)
        self._u = u
    
    

class SuctionCP:
    
    def __init__(self, uDist):
             
                        
        
                    
        