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
     """
     Generates the 3D radial turbine from a given 2D meridional channel and 
     information about the inflow and outflow angles. 
     """
    
    def __init__(self,Meridional_Channel_2D,beta_dist_hub=[20.0,15.0,-10.0,-50.0],beta_dist_tip=[20.0,5.0,-20,-51.0]):
        self._length = Meridional_Channel_2D._length
        self._beta_dist_hub = beta_dist_hub
        self._beta_dist_tip = beta_dist_tip
        self._cps_low_2D = Meridional_Channel_2D._cps_low
        self._cps_up_2D  = Meridional_Channel_2D._cps_up
        self._hub_2D = Bspline(self._cps_low_2D)
        self._tip_2D = Bspline(self._cps_up_2D)
        self._hub_3D,self._radius_hub = self._3D_curve(self._hub_2D,beta_dist_hub)
        self._tip_3D,self._radius_tip = self._3D_curve(self._tip_2D,beta_dist_tip)

    def _3D_curve(self,curve,beta_dist):
        length = self._length
        beta_interpolated = Interpolation_linear(beta_dist)()
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

        
class Interpolation_linear:
    """
    The class accepts a list of values, finds 'n' interpolation values between
    the original values and returns the list. 
    Total added values = n * (amount_of_values - 1)
    """    
    
    def __init__(self,values,n = 10):
        self._values = values
        self._n = n # number of interpolations between values   
        
    def __call__(self):    
        self._values_interpolated = self._interpolation()
        return self._values_interpolated
        
    def _interpolation(self):
        values, n = self._values, self._n
        
        amount = len(values)
        values_interpolated = zeros(amount+amount*n-n)
        
        i = 0 
        while i < amount-1:
            ii = i*(n+1)
            for j in range(0,n+2):
                values_interpolated[ii+j] = (values[i]*(n+1.0-j) + values[i+1]*j)/(n+1.0)
            i += 1
        return values_interpolated
        
    def get_values(self):
        return self._values
        

if __name__=='__main__':
        
    close('all')    

    test2Dmc = MeridionalChannel(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it)  
    test = RadialTurbine(test2Dmc,beta_dist_hub,beta_dist_tip)
    testhub = Bspline(test._hub_3D)
    testtip = Bspline(test._tip_3D)
#    checkz = RadialTurbine(inlet_diameter,outlet_mean_diameter,length,R1,R2,outlet_angle,it,[0.0,0.0,0.0,-5.0],beta_dist_tip) 
#    checkzz = Bspline(checkz._hub_3D)    
    
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
    
