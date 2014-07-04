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
from needle import NeedleDef, Tip, Shaft, Corner
import preproc_prof as pp



if __name__ == '__main__':
    filename = 'axial_nozzle/div_geometry/my_data.txt'           
    div_part = pp.read_xy_p(filename)
#     print div_part
    plot(div_part[0,:], div_part[1,:], '-g')
    #this P comes from the fitting
    P = [Point(0.0, 1.0), Point(1.86507936508, 1.0), Point(3.73015873016, 1.36696824976), Point(5.59523809524, 2.79153153735), Point(7.46031746032, 2.43521509083), Point(9.3253968254, 2.805522906), Point(11.1904761905, 2.86509988565)]
    
    div2 = Bspline(P, p = 6)
#     div2.plot()
    
    fl = 0.5 
    lc = fl*div_part[0,len(div_part[0,:])-1]
    fA = 0.9
    Ac = fA*div_part[1,len(div_part[0,:])-1]
    convDef = ConvDef(l = lc , A = Ac, m = 0.0, ncp = 4)
    conv = NozzleComp(convDef)
    conv.plot()
    
    needleDef = NeedleDef(tip_pos = 5.0, del_D = 0.5, alpha = 30.0, r_curv = 0.2, r_curv_c = 0.1, length = 6.2 )
    tip  = Tip(needleDef)
    tip.plot()
    corner  = Corner(needleDef)
    corner.plot()
    
    shaft  = Shaft(needleDef)
    shaft.plot()
    
    
##   printing to the file the geometry  
    u = np.linspace(0.0, 1.0, 50)
    t = np.linspace(0.0, 1.0, 5)
    out_file = open("axial_nozzle/convdiv_geometry/nozzle.dat","w")
#     out_file.write("profile \n")
#     out_file.write("\n")
        
    for  i in xrange(50):
        x,y,z = conv(u[i])
        out_file.write("%f %f %f \n" %(x,y, 0.000))
    for  i in xrange(len(div_part[0,:])-1):
        out_file.write("%f %f %f \n" %(div_part[0,i+1],div_part[1,i+1], 0.000))
     
    out_file.close() 
    
    out_file = open("axial_nozzle/convdiv_geometry/needle.dat","w")
#     out_file.write("profile \n")
#     out_file.write("\n")
        
    for  i in xrange(50):
        x,y,z = tip(u[i])
        out_file.write("%f %f %f \n" %(x,y, 0.000))
    for  i in xrange(5):
        x,y,z = corner(t[i])
        out_file.write("%f %f %f \n" %(x,y, 0.000))    
    out_file.write("%f %f %f \n" %(div_part[0,len(div_part[0,:])-1], y, 0.000))
     
    out_file.close() 
    
       
      
#     div.plot()
#     print div.getCP()
#     
#     testP = [Point(0.0, 1.0, 0.0, 1.0), Point(0.391547869639, 1.0, 0.0, 1.0), Point(1.527864045, 2.3, 0.0, 1.0), Point(3.29771798166, 2.5, 0.0, 1.0), Point(5.527864045, 2.7527864045, 0.0, 1.0), Point(8.0, 3.0, 0.0, 1.0)]


    axis('equal')
    show() 