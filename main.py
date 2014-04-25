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
from bspline import Bspline, DerBspline, BsplineSurface
from point import Point
from profile1 import Profile1dDist, Profile1CP
# from mpl_toolkits.mplot3d import Axes3D
    


if __name__ == '__main__':

    
    camb_def = CamberlineDef(beta_in = 30.0, beta_out = -50.0, stagger = 20 )
    camb_def2= CamberlineDef()
    camb_cp  = CamberlineCP(camb_def, opt = 1)
#     camb_cp_pup  = CamberlineCP(camb_def, opt = 1, pitch = 0.5)
#     camb_cp_pdown  = CamberlineCP(camb_def, opt = 1, pitch = -0.5)
    camb_cp_h = CamberlineCP(camb_def2, opt = 1, height = 1.0)
#     pitch = 1.0
#     Ptr = Point(0.0, pitch)
    camb = Camberline(camb_cp, p = 2)
#     camb_pup = Camberline(camb_cp_pup, p = 2)
#     camb_pdown = Camberline(camb_cp_pdown, p = 2)
    camb_h = Camberline(camb_cp_h, p = 2)
    #plot a 2D periodi domain
#     camb.plotCamb()
#     camb_pup.plotCamb()
#     camb_pdown.plotCamb()
    


#     camb.getDerCamb().plot()
    udist = CamberlineUdist(5)
#     check = Bspline(camb_cp(), p = 3)
#     check_der = DerBspline(check, kth = 1)
#     check_der.plot()
    cambNormalP = CamberlineNormalP(camb, udist)
    cambNormalP_h = CamberlineNormalP(camb_h, udist)
#     print camb_cp()
#     print cambNormalP.getCambP()
#     print cambNormalP.getDerCambP()
# #     print cambNormalP.getNormalCambP()
#     cambNormalP.plotLine()
#     print udist()
     
    s_dDistCP = [Point(0.0, 0.07), Point(0.25, 0.3), Point(0.5, 0.02), Point(1.0, 0.02)]
    p_dDistCP = [Point(0.0, 0.07), Point(0.10, 0.15),Point(0.5, 0.0), Point(1.0, 0.02)]
    s_dDist = Profile1dDist(s_dDistCP)
    p_dDist = Profile1dDist(p_dDistCP)
    s_CP = Profile1CP(cambNormalP, s_dDist)
    p_CP = Profile1CP(cambNormalP, p_dDist, sp = -1.0 )
    s_CP_h = Profile1CP(cambNormalP_h, s_dDist)
    p_CP_h = Profile1CP(cambNormalP_h, p_dDist, sp = -1.0 )
#     
    s_prof = Bspline(s_CP.getProfCP(), p = 3)
    p_prof = Bspline(p_CP.getProfCP(), p = 3)
    s_prof_h = Bspline(s_CP_h.getProfCP(), p = 3)
    p_prof_h = Bspline(p_CP_h.getProfCP(), p = 3)
    s_surface = BsplineSurface([s_CP.getProfCP(), s_CP_h.getProfCP()])
    p_surface = BsplineSurface([p_CP.getProfCP(), p_CP_h.getProfCP()])
#     P1 = []
#     P2 = []
#     for i in xrange(len(s_CP.getProfCP())):
#         p = s_CP.getProfCP()[i]
#         t = p_CP.getProfCP()[i]
#         P1.append(p+Ptr)
#         P2.append(t+Ptr)
#     s_prof_tr = Bspline(P1, p = 4)
#     p_prof_tr = Bspline(P2, p = 4)    
#     s_prof.plot()
#     p_prof.plot()
#     s_prof_tr.plot()
#     p_prof_tr.plot()



#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     u = np.linspace(0.0, 1.0, 100)
# #     v = np.linspace(0.0, 1.0, 100)
#     
#     x1 = np.ones((100), dtype = float)
#     y1 = np.ones((100), dtype = float)
#     z1 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x1[i], y1[i], z1[i] = camb.getCamb(u[i])
#     
#     x2 = np.ones((100), dtype = float)
#     y2 = np.ones((100), dtype = float)
#     z2 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x2[i], y2[i], z2[i] = camb_h.getCamb(u[i])
#     
#     x3 = np.ones((100), dtype = float)
#     y3 = np.ones((100), dtype = float)
#     z3 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x3[i], y3[i], z3[i] = s_prof(u[i])
#      
#      
#     x4 = np.ones((100), dtype = float)
#     y4 = np.ones((100), dtype = float)
#     z4 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x4[i], y4[i], z4[i] = p_prof(u[i])
#                  
#     x5 = np.ones((100), dtype = float)
#     y5 = np.ones((100), dtype = float)
#     z5 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x5[i], y5[i], z5[i] = s_prof_h(u[i])
#     x6 = np.ones((100), dtype = float)
#     y6 = np.ones((100), dtype = float)
#     z6 = np.ones((100), dtype = float)
#     for i in xrange(100):
#         x6[i], y6[i], z6[i] = p_prof_h( u[i])        
#     
#     ax.plot(x1, y1, z1, '-b')
#     ax.plot(x2, y2, z2, '-b')
#     ax.plot(x3, y3, z3, '-b')
#     ax.plot(x4, y4, z4, '-b')
#     ax.plot(x5, y5, z5, '-b')
#     ax.plot(x6, y6, z6, '-b')
#     ax.scatter(prova.get_x()[1][:], prova.get_y()[1][:], prova.get_z()[1][:])
#     ax.scatter(prova.get_x()[0][:], prova.get_y()[0][:], prova.get_z()[0][:])
#     ax.plot(prova.get_x()[1][:], prova.get_y()[1][:], prova.get_z()[1][:], '-r')
#     ax.plot(prova.get_x()[0][:], prova.get_y()[0][:], prova.get_z()[0][:], '-r')



# plot a 3D profile    
# 
#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     u = np.linspace(0.0, 1.0, 100)
#     v = np.linspace(0.0, 1.0, 100)
#      
#     xs = np.ones((100, 100), dtype = float)
#     ys = np.ones((100, 100), dtype = float)
#     zs = np.ones((100, 100), dtype = float)
#     xp = np.ones((100, 100), dtype = float)
#     yp = np.ones((100, 100), dtype = float)
#     zp = np.ones((100, 100), dtype = float)
#     for i in xrange(100):
#         for j in xrange(100):
#             xs[i][j], ys[i][j], zs[i][j] = s_surface(v[i],u[j])
#             xp[i][j], yp[i][j], zp[i][j] = p_surface(v[i],u[j])    
#             
#         ax.plot(xs[i][:], ys[i][:], zs[i][:], '-b')
#         ax.plot(xp[i][:], yp[i][:], zp[i][:], '-b')
#         
#     ax.scatter(s_surface.get_x()[1][:], s_surface.get_y()[1][:], s_surface.get_z()[1][:])
#     ax.scatter(s_surface.get_x()[0][:], s_surface.get_y()[0][:], s_surface.get_z()[0][:])
#     ax.plot(s_surface.get_x()[1][:], s_surface.get_y()[1][:], s_surface.get_z()[1][:], '-r')
#     ax.plot(s_surface.get_x()[0][:], s_surface.get_y()[0][:], s_surface.get_z()[0][:], '-r')
#     ax.scatter(p_surface.get_x()[1][:], p_surface.get_y()[1][:], p_surface.get_z()[1][:])
#     ax.scatter(p_surface.get_x()[0][:], p_surface.get_y()[0][:], p_surface.get_z()[0][:])
#     ax.plot(p_surface.get_x()[1][:], p_surface.get_y()[1][:], p_surface.get_z()[1][:], '-r')
#     ax.plot(p_surface.get_x()[0][:], p_surface.get_y()[0][:], p_surface.get_z()[0][:], '-r')     

    
    axis('equal')
    show()
    
    
    