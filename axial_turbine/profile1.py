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
import camberline
from camberline import *



    

class Profile1dDist:
    
    def __init__(self, dDistCP, p = 1):
        self._dDistCP = dDistCP
        self._p = p
        self._dDist = Part.BSplineCurve()
#         print dDistCP
        self._dDist.buildFromPoles(dDistCP, False, p)
    
    def __call__(self):
        return self._dDist
           
        
class Profile1CP:
    
    def __init__(self, cambTangentP, dDists, dDistp, p = 3):              
        self._cambTangentP = cambTangentP
        self._dDists = dDists
        self._dDistp = dDistp
        self._p = p
        self._calcdDistP()
        self._calcProfCP()
        
        
    def __call__(self):
        return 0.0
    
    
    def _calcdDistP(self):
        dDists, dDistp,  uDist = self._dDists,self._dDistp, self._cambTangentP.getCambUdist()
        dDistPs = []
        dDistPp = []
        for i in xrange(len(uDist)):
            dDistPs.append(dDists().value(uDist[i]))
            dDistPp.append(dDistp().value(uDist[i]))
        self._dDistPs = dDistPs
        self._dDistPp = dDistPp
                
        
        
    def _calcProfCP(self):
        chord, cambTangentP, dDistPs, dDistPp = self._cambTangentP.getCamb().getCambDef().getC(), self._cambTangentP, self._dDistPs, self._dDistPp
        profCPs = []
        profCPp = []
        normT = cambTangentP.getTangP()[0][0].normalize()
#         print cambTangentP.getTangP()[0][0].normalize()
        normalPp = Base.Vector(-normT.y, normT.x, 0.0)
        normalPs = normalPp.negative()
        profCPs.append(cambTangentP.getCambP()[0])
        profCPp.append(cambTangentP.getCambP()[0])
        print len(dDistPs)
        for i in xrange(len(dDistPs)):
            normT = cambTangentP.getTangP()[0][0].normalize()
#         print cambTangentP.getTangP()[0][0].normalize()
            P = cambTangentP.getCambP()[i]
            normT = cambTangentP.getTangP()[i][0].normalize()
            normalPs = Base.Vector(-normT.y, normT.x, 0.0)
            normalPp = normalPs.negative()
            profCPs.append(P.add(normalPs.multiply(dDistPs[i].y*chord)))
            profCPp.append(P.add(normalPp.multiply(dDistPp[i].y*chord)))
        profCPs.append(cambTangentP.getCambP()[len(dDistPs)-1])
        profCPp.append(cambTangentP.getCambP()[len(dDistPs)-1])    
        self._profCPs = profCPs 
        self._profCPp = profCPp
        
    def getProfCPs(self):
        return self._profCPs
    def getProfCPp(self):
        return self._profCPp              
#         
         
        
        
#         
class Profile1:
      
    def __init__(self, cp,  p = 3, pitch = 0.0, height = 0.0 ):
        self._p = p
        self._pitch = pitch
        self._height = height
#         if height == 0.0 and pitch == 0.0:
        self._cp_s = cp.getProfCPs()
        self._cp_p = cp.getProfCPp()    
#         else:
#             self._cp_s = self._redefCP(cp_s)
#             self._cp_p = self._redefCP(cp_p)
         
        self._prof_s = Part.BSplineCurve()    
        self._prof_p = Part.BSplineCurve()
        self._prof_s.buildFromPoles(self._cp_s, False, p)
        self._prof_p.buildFromPoles(self._cp_p, False, p)
             
         
#     def _redefCP(self, cp):
#         pitch, height = self._pitch, self._height
#         cp_new = []
#         for i in xrange(len(cp)):
#             cp_new.append(Point(cp[i].x, cp[i].y + pitch, cp[i].z + height))
#         return cp_new    
         
         
    def getProfS(self):
        return self._prof_s              
     
    def getProfP(self):
        return self._prof_p
     

class SimulationDomain(object):
    
    def __init__(self, cambDef, camb, pitch = 1.0 ):
        self.cambDef = cambDef
        self.camb = camb
        self.pitch = pitch
        self._inout()
        print self.camb.getCamb().getPoles()
        self.cambPerTop = Part.BSplineCurve()
        self.cambPerTop.buildFromPoles(self.camb.getCamb().getPoles(), False, self.camb.getCamb().Degree)
        self.cambPerBot = Part.BSplineCurve()
        self.cambPerBot.buildFromPoles(self.camb.getCamb().getPoles(), False, self.camb.getCamb().Degree)
        self.cambPerTop.translate(Base.Vector(0.0, pitch/2.0, 0.0))
        self.cambPerBot.translate(Base.Vector(0.0, -pitch/2.0, 0.0))
#         print self.cambPerBot.getPoles()
#         print self.cambPerTop.getPoles()
        self.channel = [self.inlet.toShape(), self.inPerTop.toShape(), self.cambPerTop.toShape(), 
                        self.outPerTop.toShape(), self.outlet.toShape(), self.outPerBot.toShape(), 
                        self.cambPerBot.toShape(), self.inPerBot.toShape()]
       
        
    def _inout(self):
        cambDef, pitch = self.cambDef, self.pitch
        le = cambDef.getLe()
        te = cambDef.getTe()
        c_ax = cambDef.getC_ax()
        self.inlet = Part.Line(Base.Vector(le.x - c_ax, le.y - pitch/2.0, 0.0 ) , Base.Vector(le.x - c_ax, le.y + pitch/2.0, 0.0 ))
        self.outlet = Part.Line(Base.Vector(te.x + c_ax, te.y + pitch/2.0, 0.0 ) , Base.Vector(te.x + c_ax, te.y - pitch/2.0, 0.0 )) 
        self.inPerTop = Part.Line(Base.Vector(le.x - c_ax, le.y + pitch/2.0, 0.0 ) , Base.Vector(le.x, le.y + pitch/2.0, 0.0 ))
        self.inPerBot = Part.Line(Base.Vector(le.x, le.y - pitch/2.0, 0.0 ) , Base.Vector(le.x - c_ax, le.y - pitch/2.0, 0.0 ))   
        self.outPerTop = Part.Line(Base.Vector(te.x, te.y + pitch/2.0, 0.0 ) , Base.Vector(te.x + c_ax, te.y + pitch/2.0, 0.0 ))
        self.outPerBot = Part.Line(Base.Vector(te.x + c_ax , te.y - pitch/2.0, 0.0 ) , Base.Vector(te.x, te.y - pitch/2.0, 0.0 ))

if __name__=='__main__':
    
    
    cambdef= CamberlineDef(beta_in=30.0, stagger=30.0)           
    cambCP = CamberlineCP(cambdef)
    camb = Camberline(cambCP)
    udist = CamberlineUdist(5)
    cambTangentP = CamberlineTangentP(camb, udist)
    s_dDistCP = [Base.Vector(0.0, 0.07, 0.0),  Base.Vector(0.50, 0.15, 0.0),  Base.Vector(0.7, 0.02, 0.0),  Base.Vector(1.0, 0.02, 0.0)]
    p_dDistCP = [Base.Vector(0.0, 0.07, 0.0), Base.Vector(0.05, 0.10, 0.0), Base.Vector(0.5, 0.0, 0.0), Base.Vector(1.0, 0.02, 0.0)]
    s_dDist = Profile1dDist(s_dDistCP)
    p_dDist = Profile1dDist(p_dDistCP)
    profCP = Profile1CP(cambTangentP,s_dDist, p_dDist)
    profile = Profile1(profCP)
#     b2bsec = Part.makeRuledSurface(profile.getProfS().toShape(), profile.getProfP().toShape())
    b2bsec = Part.makeFilledFace([profile.getProfS().toShape(), profile.getProfP().toShape()])
    blade = b2bsec.extrude(Base.Vector(0.0,0.0,1.0))
    simDom = SimulationDomain(cambdef, camb)
    dom = Part.makeFilledFace(simDom.channel)
    simdom2D = dom.cut(b2bsec)
    simdom3D = simdom2D.extrude(Base.Vector(0.0,0.0,1.0))
    
    
#     edge1 = Part.Edge(profile.getProfS())
#     edge2 = Part.Edge(profile.getProfP())
#     wprofile = Part.Wire(profile.getProfS().toShape().Edges[0],profile.getProfP().toShape().Edges[0])
#     a = profile.getProfP().toShape()


    

    FreeCADGui.showMainWindow()
#     Part.show(wprofile)
#     Part.show(profile.getProfS().toShape())
#     Part.show(profile.getProfP().toShape())
    Part.show(blade)
#     Part.show(b2bsec)
#     Part.show(simDom.inlet.toShape())
#     Part.show(simDom.outlet.toShape())
#     Part.show(simDom.inPerBot.toShape())
#     Part.show(simDom.inPerTop.toShape())
#     Part.show(simDom.outPerTop.toShape())
#     Part.show(simDom.outPerBot.toShape())
#     Part.show(simDom.cambPerBot.toShape())
#     Part.show(simDom.cambPerTop.toShape())
#     Part.show(camb.getCamb().toShape())
    Part.show(simdom3D)  
        
#     Part.show(camb._cambCP._teLine.toShape())
#     Part.show(camb._cambCP._leLine.toShape())
    FreeCADGui.exec_loop()            
        
             
                        
        
                    
        