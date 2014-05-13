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


    convDef = ConvDef()
    divDef = DivDef()
    conv = ConvDiv(convDef)
    div = ConvDiv(divDef)
    
    conv.plot()
    div.plot()

    axis('equal')
    show() 