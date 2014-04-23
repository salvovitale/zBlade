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
from point import Point, Angle, normVector, normal2D, Vector
from bspline import Bspline
from line import Line
from math import sqrt, pi, cos, sin, tan


    

class Profile1dDist:
    
    def __init__(self, dDistCP, p = 2):
        self._dDistCP = dDistCP
        self._p = p
        self._dDist = Bspline(dDistCP,p)
    
    def __call__(self, u):
        return self._dDist(u)
        
    
    def plot(self):
        self._dist.plot()   
        
class Profile1CP:
    
    def __init__(self, camberNormalP, dDist, p = 3, sp = 1.0):              
        self._camberNormalP = camberNormalP
        self._dDist = dDist
        self._p = p
        self._sp = sp
        self._calcdDistP()
        self._calcProfCP()
        
        
    def __call__(self):
        return 0.0
    
    
    def _calcdDistP(self):
        dDist, uDist = self._dDist, self._camberNormalP.getCambUdist()
        dDistP = np.ones(len(uDist))
        for i in xrange(len(uDist)):
            x, dDistP[i], z = dDist(uDist[i])
        self._dDistP = dDistP
                
        
        
    def _calcProfCP(self):
        camberNormalP, dDistP, sp = self._camberNormalP, self._dDistP, self._sp
        profCP = []
        profCP.append(camberNormalP.getNormalLine()[0](0.0))
        for i in xrange(len(dDistP)):
            profCP.append(camberNormalP.getNormalLine()[i](dDistP[i]*sp))
        profCP.append(camberNormalP.getNormalLine()[len(dDistP)-1](0.0)) 
        self._profCP = profCP
        
    def getProfCP(self):
        return self._profCP              
        
         
        
        
        
        
            
        
        

        
            
        
             
                        
        
                    
        