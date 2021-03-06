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
        chord, camberNormalP, dDistP, sp = self._camberNormalP.getCamb().getCambDef().getC(), self._camberNormalP, self._dDistP, self._sp
        profCP = []
        profCP.append(camberNormalP.getNormalLine()[0](0.0))
        for i in xrange(len(dDistP)):
            profCP.append(camberNormalP.getNormalLine()[i](dDistP[i]*sp*chord))
        profCP.append(camberNormalP.getNormalLine()[len(dDistP)-1](0.0)) 
        self._profCP = profCP
        
    def getProfCP(self):
        return self._profCP              
        
         
        
        
        
class Profile1:
     
    def __init__(self, cp_s, cp_p, p = 3, pitch = 0.0, height = 0.0 ):
        self._p = p
        self._pitch = pitch
        self._height = height
        if height == 0.0 and pitch == 0.0:
            self._cp_s = cp_s
            self._cp_p = cp_p
        else:
            self._cp_s = self._redefCP(cp_s)
            self._cp_p = self._redefCP(cp_p)
        
        self._prof_s = Bspline(self._cp_s, p = p)    
        self._prof_p = Bspline(self._cp_p, p = p)
            
        
    def _redefCP(self, cp):
        pitch, height = self._pitch, self._height
        cp_new = []
        for i in xrange(len(cp)):
            cp_new.append(Point(cp[i].x, cp[i].y + pitch, cp[i].z + height))
        return cp_new    
        
        
    def getProfS(self, u):
        return self._prof_s(u)              
    
    def getProfP(self, u):
        return self._prof_p(u)
    
    def plot(self):
        self._prof_s.plot()        
        self._prof_p.plot()
        

        
            
        
             
                        
        
                    
        