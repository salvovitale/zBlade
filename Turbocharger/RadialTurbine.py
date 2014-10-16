# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 11:53:10 2014

@author: Roel
"""

zblade_folder= '../'
sys.path.insert(0, zblade_folder)

from bspline import *
from Meridional import *
from math import pi, sin, cos, tan, asin
from point import *
#import numpy as np

if __name__=='__main__':
#==============================================================================
    # Input
    inlet_diameter = 1.0 # meridional channel inlet_diameter
    outlet_mean_diameter = 0.5
    length = 0.7 # meridional channel length
    R1 = 0.2 # inlet width meridional channel
    R2 = 0.4 # outlet width meridional channel
    outlet_angle = 0 # channel outlet angle in degrees
    beta_dist_hub = [20.0,5.0,-20.0,-47.0] 
    beta_dist_tip = [20.0,5.0,-30,-51.0]
#    beta_dist_hub = [2.0,8.0,16.0,36.0] #
    # thickness distribution (constant?)
    it = 1 # number of iterations
#==============================================================================  
    
class RadialTurbine:
    
    def __init__(self,inlet_diameter=1,outlet_mean_diameter=0.5,length=0.7,R1=0.2,R2=0.4,
                 outlet_angle=0,it=1,beta_dist_hub=[20.0,15.0,-10.0,-50.0],beta_dist_tip=[20.0,5.0,-20,-51.0]):
        self._inlet_diameter = inlet_diameter
        self._length = length
        self._R1 = R1
        self._R2 = R2
        self._outlet_angle = outlet_angle
        self._it = it
        self._beta_dist_hub = beta_dist_hub
        self._beta_dist_tip = beta_dist_tip
        self._cps_low_2D = MeridionalChannel(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it)._cps_low
        self._cps_up_2D = MeridionalChannel(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it)._cps_up
        self._hub_2D = Bspline(self._cps_low_2D)
        self._tip_2D = Bspline(self._cps_up_2D)
        self._hub_3D,self._radius_hub = self._3D_curve(self._hub_2D,beta_dist_hub)
        self._tip_3D,self._radius_tip = self._3D_curve(self._tip_2D,beta_dist_tip)
  

    def _3D_curve(self,curve,beta_dist):
        length = self._length
        beta_interpolated = self._beta_interpolate(beta_dist)
        nbeta = len(beta_interpolated)
        x = zeros(nbeta)
        y = ones(nbeta)
        z = zeros(nbeta)
        radius = zeros(nbeta)
        theta = zeros(nbeta)
#        rcheck = ones(nbeta)
        x[0] = curve(0)[0]
        y[0] = curve(0)[1]
        z[0] = curve(0)[2]
        radius[0] = curve(0)[1]
        # find 3D x,y,z coordinates 
        for i in range(1,nbeta):
            x[i] = curve(i/(nbeta-1.0))[0]
            radius[i] = curve(i/(nbeta-1.0))[1]
            theta[i] = tan((beta_interpolated[i]-beta_interpolated[i-1])*pi/180.0)*(x[i]-x[i-1])/radius[i]+theta[i-1]
##            print theta*180.0/pi
            z[i] = radius[i]*sin(theta[i])+z[i-1] #-theta[i-1]
#            z[i] = (i/(nbeta-1.0))-(i-1)/(nbeta-1.0)*tan(Angle(beta_interpolated[i])())+z[i-1]
#            z[i] = (x[i]-x[i-1])*tan(Angle(beta_interpolated[i])())+z[i-1]
            y[i] = sqrt(radius[i]**2-z[i]**2)
#            rcheck[i] = sqrt(y[i]**2+z[i]**2)
#        print 'radius = %s '%radius
        print 'theta = %s' %(theta*180.0/pi)
#        print 'x = %s' %x
        print 'y = %s' %y
#        print 'z = %s' %z
#        print 'rcheck = %s' %rcheck
        cps_3D = list(zeros(len(x)))
        for i in range(0,len(cps_3D)):
            cps_3D[i] = Point(x[i],y[i],z[i])
        print cps_3D
            
        return cps_3D,radius

    def _beta_interpolate(self,beta_dist,n=12):
        # n = number of (linear) interpolations between betas
        nbeta = len(beta_dist)
        beta_dist_interpolated = zeros(nbeta+nbeta*n-n)
        
        i = 0 
        while i < nbeta-1:
#            print i
            ii = i*(n+1)
#            print ii
            for j in range(0,n+2):
                beta_dist_interpolated[ii+j] = (beta_dist[i]*(n+1-j) + beta_dist[i+1]*j)/(n+1)
#            print beta_dist_interpolated
            i += 1
#        print beta_dist_interpolated
        return beta_dist_interpolated

if __name__=='__main__':
        
    close('all')    
    test = RadialTurbine(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it,beta_dist_hub,beta_dist_tip)
    testhub = Bspline(test._hub_3D)
    testtip = Bspline(test._tip_3D)
    checkz = RadialTurbine(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it,[0.0,0.0,0.0,-5.0],beta_dist_tip) 
    checkzz = Bspline(checkz._hub_3D)    
    
    from mpl_toolkits.mplot3d import Axes3D
    fig1 = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlim3d(0,1)
    ax.set_ylim3d(0,1)
    ax.set_zlim3d(-0.5,0.5)
    xlabel('x')
    ylabel('y')
    ax.plot(testhub._X,testhub._Y,testhub._Z,'b')
    ax.plot(testtip._X,testtip._Y,testtip._Z,'g')
    
    fig2 = plt.figure()
    plt.plot(testhub._X,testhub._Y,'b')
    plt.plot(testtip._X,testtip._Y,'g')
    plt.plot(testhub._X,test._radius_hub,'r')
    xlabel('x')
    ylabel('y')
    
    fig3 = plt.figure()
    plt.plot(testhub._Z,testhub._Y)
    plt.plot(testtip._Z,testtip._Y)
#    plt.plot(checkzz._Z,checkzz._Y)
    axis([-1,0.1,-0.1,1])
    xlabel('z')
    ylabel('y')
#    testhub.plot()
#    testtip.plot()
    
    #pakt z niet goed mee met plotten?
    #interpoleren tussen hub en tip
    #andere tip beta angles
    
#    for i in range(1,nbeta):
#            x[i] = (i/(nbeta-1.0))*(curve(1)[0]-curve(0)[0])+curve(0)[0]
#            radius[i] = curve(x[i]/(curve(1)[0]-curve(0)[0]+curve(0)[0]))     [1]
##            theta[i] = tan((beta_interpolated[i]-beta_interpolated[i-1])*pi/180.0)*(x[i]-x[i-1])/radius[i]+theta[i-1]
##            print theta*180.0/pi
##            z[i] = radius[i]*sin(theta[i]-theta[i-1])+z[i-1]
#            z[i] = (x[i]-x[i-1])*tan(Angle(beta_interpolated[i])())+z[i-1]
#            y[i] = sqrt(radius[i]**2-z[i]**2)
##            rcheck[i] = sqrt(y[i]**2+z[i]**2)