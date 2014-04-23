"""
Created on Apr 19, 2014

@author: salvatorevitale
"""
import sys
zblade_folder= '../'
sys.path.insert(0, zblade_folder)
from matplotlib.pylab import *
# from point import Point, Norm, Angle
# from line  import Line
from point import  Point, Angle, normVector, normal2D, Vector
from bspline import *
from line import Line
from bspline import Bspline, DerBspline
from math import sqrt, pi, cos, sin, tan

  


class CamberlineDef:
    
    def __init__(self, le = 0.0, beta_in = 30.0 , beta_out = -70.0 , c_ax = 1.0, stagger = 40.0):
        if le == 0.0:
            self._le = Point(0.0, 0.0)
        else:
            self._le = le
        self._beta_in  = Angle(beta_in)
        self._beta_out = Angle(beta_out)
        self._c_ax = c_ax
        self._stagger = Angle(stagger)
        self._calc_te()
    
    def _calc_te(self):
        le, c_ax, stagger = self._le, self._c_ax, self._stagger
        self._te = Point(le.x + c_ax, le.y - c_ax*tan(stagger()))
    
    def getLe(self):
        return self._le    
    
    def getBin(self):
        return self._beta_in()
    
    def getBout(self):
        return self._beta_out()
    
    def getC_ax(self):
        return self._c_ax
    
    def getStagger(self):
        return self._stagger
    
    def getTe(self):
        return self._te
    

class CamberlineCP:
    
    def __init__(self, cambDef, opt = 1, t1 = 0.5, t2 = 0.5):
        self._cambDef = cambDef
        self._opt = opt
        self._t1 = t1
        self._t2 = t2
        #definire le due linee per calcolare l'intersezione
        self._leLine = Line(cambDef.getLe(), Vector(cos(cambDef.getBin()),sin(cambDef.getBin())))
        self._teLine = Line(cambDef.getTe(), Vector(cos(cambDef.getBout()),sin(cambDef.getBout())))
        self._interPoint()
        self._controlPoint()
    
    def __call__(self):
        return self._cp
        
    def _interPoint(self):
        self._int_p = self._leLine.intersect(self._teLine)
              
    def _controlPoint(self):
        le, te, int_p, opt, t1, t2 = self._cambDef.getLe(), self._cambDef.getTe(), self._int_p, self._opt, self._t1, self._t2  
        if opt == 1:
            self._cp =[le, int_p, te]
        if opt == 2:
            self._cp = [le, le.towards(int_p, t1), te.towards(int_p, t2), te]
                
                
           
               
            
        
class Camberline:
    
    def __init__(self, cambCP, p = 2):
        self._cambCP = cambCP
        self._p = p
        self._camb = Bspline(cambCP(), p = p ) 
        self._derCamb = DerBspline(self._camb,  kth = 1)
    
    def getCamb(self, u):
        return self._camb(u)
    
    def getDerCamb(self, u):
        return self._derCamb(u)
    
    def plotCamb(self):
        self._camb.plot()
            
    def plotDerCamb(self):
        self._deramb.plot()    

class CamberlineUdist(object):
    """
    classdocs
    """


    def __init__(self, nintP, opt = 1, a = 0.0 , b = 1.0):
        
        """
        Constructor
        """
        self._opt = opt
        self._nintP = nintP
        self._a = a
        self._b = b
        self._udistribution()
        
        
    
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
            u[i] = a + d*i
        self._u = u
        
        
        
        
class CamberlineNormalP:
    
    def __init__(self, camb,  uDist):
        
        self._camb = camb
        self._uDist = uDist        
        self._calcCambP()
        self._calcDerCambP()
        self._calcNormalCambP()
        self._calcNormalLine()
       
    def __call__(self):
        
        return 0.0   
        
    def _calcCambP(self):
        camb, uDist = self._camb, self._uDist
        P = []    
        for i in xrange(len(uDist())):
            x, y, z = camb.getCamb(uDist()[i])
            P.append(Point(x, y, z))
        self._cambP = P
    
    
    def _calcDerCambP(self):
        camb, uDist = self._camb, self._uDist
        P = []    
        for i in xrange(len(uDist())):
            x, y, z = camb.getDerCamb(uDist()[i])
            P.append(Vector(x,y,z))
        self._derCambP = P
    
    def _calcNormalCambP(self):
        derCambP = self._derCambP
        normalCambP = []
        for i in xrange(len(derCambP)):
            normalCambP.append((normal2D(derCambP[i].norm())))
        self._normalCambP = normalCambP
        
    def _calcNormalLine(self):
        normalCambP, cambP = self._normalCambP, self._cambP
        L = []
        for i in xrange(len(cambP)):
            L.append(Line(cambP[i], normalCambP[i]))
        self._normalLine = L                      
        
    def getCambP(self):
        return self._cambP
    
    def getDerCambP(self):
        return self._derCambP
    
    def getNormalCambP(self):
        return self._normalCambP
    
    def getNormalLine(self):
        return self._normalLine   
    
    def getCambUdist(self):
        return self._uDist()
    
    def plotLine(self):
        normalLine = self._normalLine
        for i in xrange(len(normalLine)):
            normalLine[i].plot()
                 

# class Camberline(object):
#     """
#     classdocs
#     """
# 
# 
#     def __init__(self, beta_in, beta_out, stugger):
#         
#         """
#         Constructor
#         
#         """



if __name__=='__main__':
    
#     prova = Angle(30.0)
#     print prova.getAngle()
#     print prova()
#     
    prova = CP_camberline()
    camb = Bspline(prova())
    camb.plot()
    show()

        