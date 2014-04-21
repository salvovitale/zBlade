'''
Created on Apr 20, 2014

@author: salvatorevitale
'''
import sys
# from zBlade.axial_turbine.camberline import CamberlineCP
axial= 'axial_turbine/'
sys.path.insert(0, axial)
from matplotlib.pylab import *
from camberline import CamberlineCP, CamberlineDef 
from bspline import Bspline


if __name__ == '__main__':
    camb_def = CamberlineDef()
    camb_cp  = CamberlineCP(camb_def, opt = 2)
    camb = Bspline(camb_cp(), 3)
    camb.plot()
    axis('equal')
    show()
    
    