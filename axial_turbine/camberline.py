"""
Created on Apr 19, 2014

@author: salvatorevitale
"""
import sys
zblade_folder= '../'
sys.path.insert(0, zblade_folder)
# FREECADPATH = 'D:\FreeCAD_014\lib' # path to your FreeCAD.so or FreeCAD.dll file
# sys.path.append(FREECADPATH)

from matplotlib.pylab import *
# from point import Point, Norm, Angle
# from line  import Line
import FreeCAD
import FreeCADGui
# from FreeCADGui import UiLoader 
import Part
from FreeCAD import Base, Vector
import Draft
import DraftGeomUtils
import PySide
from PySide import QtCore,QtGui

# import Workbench

# from PyQt4 import QtGui
# from PyQt4 import QtCore
      
   


   




# from FreeCADBase import Vector


 
  

from point import  Angle, normVector, normal2D
# from bspline import *
# from line import Line
# from bspline import Bspline, DerBspline
from math import sqrt, pi, cos, sin, tan

  


class CamberlineDef(object):
    """
    I think is better to bring pithc and height directly in camberlineCP
    """
    
    def __init__(self, le = 0.0, beta_in = 0.0 , beta_out = -70.0 , c_ax = 1.0, stagger = 40.0):
        if le == 0.0:
            self._le = Base.Vector(0.0, 0.0)
        else:
            self._le = Base.Vector(le.x, le.y)
        self._beta_in  = Angle(beta_in)
        self._beta_out = Angle(beta_out)
        self._c_ax = c_ax
        self._stagger = Angle(stagger)
        self._calc_te()
        self._calc_chord()
         
    
    def _calc_te(self):
        le, c_ax, stagger= self._le, self._c_ax, self._stagger
        self._te = Base.Vector(le.x + c_ax, le.y - c_ax*tan(stagger()))
    
    def _calc_chord(self):
        le, te = self._le, self._te
        self._c = le.distanceToPoint(te)
        print(self._c)
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
    
    def getC(self):
        return self._c
    
    

class CamberlineCP(object):
    
    def __init__(self, cambDef, opt = 2, t1 = 0.5, t2 = 0.5,  pitch = 0.0, height = 0.0):
        self._cambDef = cambDef
        self._opt = opt
        self._t1 = t1
        self._t2 = t2
        self._pitch = pitch
        self._height = height
        self._translation()
        #definire le due linee per calcolare l'intersezione
        self._leLine = Part.Line(self._le, Base.Vector(3.0*cos(cambDef.getBin()),3.0*sin(cambDef.getBin())))
        self._teLine = Part.Line(self._te, Base.Vector(-3.0*cos(cambDef.getBout()),-3.0*sin(cambDef.getBout())))
        self._interPoint()
        
        self._controlPoint()
    
    def __call__(self):
        return self._cp
        

    def _translation(self):
        leDef, teDef, pitch, height = self._cambDef.getLe(), self._cambDef.getTe(), self._pitch, self._height
        self._le = Base.Vector(leDef.x, leDef.y + pitch, height)
        self._te = Base.Vector(teDef.x, teDef.y + pitch, height)    
    
    def _interPoint(self):
#         print self._leLine.intersect(self._teLine)
        inter_plane = Part.Plane() # (10,10, self._le)
        int_p = self._leLine.intersect2d(self._teLine, inter_plane)
        print int_p, int_p[0][0], int_p[0][1]
        self._int_p = Base.Vector(int_p[0][0], int_p[0][1])
        
     
    def _controlPoint(self):
        le, te, int_p, opt = self._le, self._te, self._int_p, self._opt  
        if opt == 1:
            self._cp =[le, int_p, te]
        if opt == 2:
            midp1 = le.add(int_p)
            midp2 = te.add(int_p)
            self._cp = [le, Base.Vector(midp1.x/2.0,midp1.y/2.0 ), Base.Vector(midp2.x/2.0,midp2.y/2.0 ), te]
#         print self._cp
    def getCambDef(self):
        return self._cambDef            
            
        
class Camberline(object):
     
    def __init__(self, cambCP, p = 2):
        self._cambCP = cambCP
        self._cambDef = cambCP.getCambDef()
        self._p = p
        self._constructCamberline()
        
    def __call__(self):
        return self._camb    
         
#         self._derCamb = DerBspline(self._camb,  kth = 1)
     
    def _constructCamberline(self):
        camb = Part.BSplineCurve()
        camb.buildFromPoles(self._cambCP(),False, self._p)
        print camb.getPoles(), camb.KnotSequence, camb.Degree
        self._camb = camb     
#         self._camb = Draft.makeBSpline(self._cambCP())
#         self._camb.Degree
#         i = 0
#         for p in self._cambCP():
#             i= i+1
#         print self._camb.getPoles()
#         print self._cambCP()[0]
#         self._camb.setPole(0, self._cambCP()[0])
#         self._camb.setPole(1, self._cambCP()[1])
#         self._camb.setPole(2, self._cambCP()[2])    
    
    def getCamb(self):
        return self._camb
     
    def getCambDef(self):
        return self._cambDef 
        

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
        
        
        
        
class CamberlineTangentP:
     
    def __init__(self, camb,  uDist):
         
        self._camb = camb
        self._uDist = uDist        
        self._calcCambP()
#         self._calcDerCambP()
        self._calcTangCambAtP()
        
        
    def __call__(self):
         
        return 0.0   
         
    def _calcCambP(self):
        camb, uDist = self._camb, self._uDist
        P = []    
        for i in xrange(len(uDist())):
            P.append(camb().value(uDist()[i]))
        self._cambP = P
        print P
     
     
#     def _calcDerCambP(self):
#         camb, uDist = self._camb, self._uDist
#         P = []    
#         for i in xrange(len(uDist())):
#             x, y, z = camb.getDerCamb(uDist()[i])
#             P.append(Vector(x,y,z))
#         self._derCambP = P
     
    def _calcTangCambAtP(self):
        camb, uDist = self._camb, self._uDist
        camb, uDist = self._camb, self._uDist
        T = []    
        for i in xrange(len(uDist())):
            T.append(camb().tangent(uDist()[i]))
        self._cambTangP = T
        print T
                            
     
    def getCamb(self):
        return self._camb
         
    def getCambP(self):
        return self._cambP
     
    def getTangP(self):
        return self._cambTangP
         
    def getCambUdist(self):
        return self._uDist()
     
             

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
#     prova = CP_camberline()
#     camb = Bspline(prova())
#     camb.plot()
#     show()
    cambdef= CamberlineDef()           
    cambCP = CamberlineCP(cambdef)
    camb = Camberline(cambCP)
    udist = CamberlineUdist(5)
    cambTangentP = CamberlineTangentP(camb, udist)

    
    
#     FreeCADGui.showMainWindow()
#     Part.show(camb._camb.toShape())
# #     Part.show(camb._cambCP._teLine.toShape())
# #     Part.show(camb._cambCP._leLine.toShape())
#     FreeCADGui.exec_loop()
    

    

                   
    
        