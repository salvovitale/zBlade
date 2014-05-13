"""
Created on May 13, 2014

@author: salvovitale
"""
import sys
import string
# from zBlade.axial_turbine.camberline import CamberlineCP
axial= 'axial_nozzle/'
sys.path.insert(0, axial)
from matplotlib.pylab import *
from nozzle2D import ConvDef, DivDef, ConvDiv
from point import Point
from bspline import Bspline



if __name__ == '__main__':


    convDef = ConvDef(l = 3.0, A = 2.0, m = -0.1, ncp = 4, type = 2)
    divDef = DivDef(l = 8.0, A = 3.0, m = 0.1, ncp = 6, type = 3)
    conv = ConvDiv(convDef)
    div = ConvDiv(divDef)
    
    conv.plot()
#     div.plot()
    print div.getCP()
    
    testP = [Point(0.0, 1.0, 0.0, 1.0), Point(0.391547869639, 1.0, 0.0, 1.0), Point(1.527864045, 2.3, 0.0, 1.0), Point(3.29771798166, 2.5, 0.0, 1.0), Point(5.527864045, 2.7527864045, 0.0, 1.0), Point(8.0, 3.0, 0.0, 1.0)]
    div2 = Bspline(testP, p = 4)
    div2.plot()

    axis('equal')
    show() 