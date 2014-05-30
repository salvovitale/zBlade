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
from nozzle2D import ConvDef, DivDef, NozzleComp
from point import Point
from bspline import Bspline
import preproc_prof as pp



if __name__ == '__main__':
    filename = 'axial_nozzle/div_geometry/my_data.txt'           
    div_part = pp.read_xy_p(filename)
    print div_part
    
    fl = 0.5 
    lc = fl*div_part[0,len(div_part[0,:])-1]
    fA = 0.9
    Ac = fA*div_part[1,len(div_part[0,:])-1]
    convDef = ConvDef(l = lc , A = Ac, m = 0.0, ncp = 4)
    conv = NozzleComp(convDef)
    conv.plot()
    plot(div_part[0,:], div_part[1,:], '-g')
    
    u = np.linspace(0.0, 1.0, 50)
    out_file = open("axial_nozzle/convdiv_geometry/test1.dat","w")
#     out_file.write("profile \n")
#     out_file.write("\n")
       
    for  i in xrange(50):
        x,y,z = conv(u[i])
        out_file.write("%f %f %f \n" %(x,y, 0.000))
    for  i in xrange(len(div_part[0,:])-1):
        out_file.write("%f %f %f \n" %(div_part[0,i+1],div_part[1,i+1], 0.000))
    
    out_file.close()    
    
#     div.plot()
#     print div.getCP()
#     
#     testP = [Point(0.0, 1.0, 0.0, 1.0), Point(0.391547869639, 1.0, 0.0, 1.0), Point(1.527864045, 2.3, 0.0, 1.0), Point(3.29771798166, 2.5, 0.0, 1.0), Point(5.527864045, 2.7527864045, 0.0, 1.0), Point(8.0, 3.0, 0.0, 1.0)]
#     div2 = Bspline(testP, p = 4)
#     div2.plot()

    axis('equal')
    show() 