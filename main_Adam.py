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



if __name__ == '__main__':


    convDef = ConvDef(l = 3.0, A = 2.0, m = 0.0, ncp = 4, type = 1)
    divDef = DivDef(l = 8.0, A = 3.0, m = 0.0, ncp = 4, type = 1)
    conv = ConvDiv(convDef)
    div = ConvDiv(divDef)
    
    conv.plot()
    div.plot()

    axis('equal')
    show() 