'''
Created on Apr 20, 2014

@author: salvatorevitale
'''
import sys

# from zBlade.axial_turbine.camberline import CamberlineCP
axial= 'axial_turbine/'
sys.path.insert(0, axial)
from matplotlib.pylab import *
from camberline import CamberlineCP, CamberlineDef, Camberline, CamberlineUdist, CamberlineNormalP  
from bspline import Bspline, DerBspline
from point import Point
from profile1 import Profile1dDist, Profile1CP


if __name__ == '__main__':
    camb_def = CamberlineDef()
    camb_cp  = CamberlineCP(camb_def, opt = 2)
    camb = Camberline(camb_cp, p = 3)
    camb.plotCamb()
#     camb.getDerCamb().plot()
    udist = CamberlineUdist(5)
#     check = Bspline(camb_cp(), p = 3)
#     check_der = DerBspline(check, kth = 1)
#     check_der.plot()
    cambNormalP = CamberlineNormalP(camb, udist)
#     print camb_cp()
#     print cambNormalP.getCambP()
#     print cambNormalP.getDerCambP()
# #     print cambNormalP.getNormalCambP()
#     cambNormalP.plotLine()
#     print udist()
    
    s_dDistCP = [Point(0.0, 0.1), Point(0.3, 0.3), Point(0.5, 0.0), Point(1.0, 0.02)]
    p_dDistCP = [Point(0.0, 0.1), Point(0.10, 0.1),Point(0.5, 0.0), Point(1.0, 0.02)]
    s_dDist = Profile1dDist(s_dDistCP)
    p_dDist = Profile1dDist(p_dDistCP)
    s_CP = Profile1CP(cambNormalP, s_dDist)
    p_CP = Profile1CP(cambNormalP, p_dDist, sp = -1.0 )
    s_prof = Bspline(s_CP.getProfCP(), p = 4)
    p_prof = Bspline(p_CP.getProfCP(), p = 4)
    s_prof.plot()
    p_prof.plot()
    axis('equal')
    show()
    
    