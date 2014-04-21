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
from point import Point, Angle, Norm
from bspline import *
from line import Line
from math import sqrt, pi, cos, sin, tan

  


class CamberlineDef:
    
    def __init__(self, le = 0.0, beta_in = 40.0 , beta_out = -60.0 , c_ax = 1.0, stagger = 30.0):
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
        self._leLine = Line(cambDef.getLe(), Norm(cos(cambDef.getBin()),sin(cambDef.getBin())))
        self._teLine = Line(cambDef.getTe(), Norm(cos(cambDef.getBout()),sin(cambDef.getBout())))
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

        