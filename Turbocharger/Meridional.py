# -*- coding: utf-8 -*-
"""
Created on Mon Oct 06 16:56:43 2014

@author: Roel
"""

zblade_folder= '../'
sys.path.insert(0, zblade_folder)

from bspline import *
from math import pi, sin, tan

if __name__=='__main__':
#==============================================================================
    # Input
    height = 1.0 # meridional channel height
    length = 1.0 # meridional channel length
    R1 = 0.2 # inlet width meridional channel
    R2 = 0.4 # outlet width meridional channel
    outlet_angle = 0 # channel outlet angle in degrees
     
    it = 1 # number of iterations
#==============================================================================           

class MeridionalChannel:
    """Constructs the 2D meridional channel and adds control points to the 
    lower and upper side of the channel's curves.
    Input:  channel height
            channel length
            channel inlet width
            channel outlet width
            channel outlet angle
            number of iterations                 
    """
    
    def __init__(self,height,length,R1,R2,outlet_angle,it):
        self._height = height
        self._length = length
        self._R1 = R1
        self._R2 = R2
        self._outlet_angle = outlet_angle
        self._it = it
        self._cps_low, self._cps_up = self.__call__()
    
    def __call__(self):
        """
        calculates start control points for upper and lower curve and activates
        _find_newcps to find the new control points. Then returns those to the 
        initiator.
        """
        height, length, R1, R2, outlet_angle, it = self._height, self._length, self._R1, self._R2, self._outlet_angle, self._it
        cpinlet = Point(0.0,height,0.0)
        cpmid = Point(0.0,-length*tan(outlet_angle*pi/180.0),0.0)
        cpoutlet = Point(length,0.0,0.0)
        cps_low = [cpinlet,cpmid,cpoutlet]
        cps_up = cps_low[:]
        cps_up[0] = Point(cps_low[0].x+R1,cps_low[0].y,   cps_low[0].z)
        cps_up[1] = Point(cps_low[1].x+R1,cps_low[1].y+R2*cos(outlet_angle*pi/180.0),cps_low[1].z)
        cps_up[2] = Point(cps_low[2].x-R2*sin(outlet_angle*pi/180.0),cps_low[2].y+R2*cos(outlet_angle*pi/180.0),cps_low[2].z)
        self._new_cps_low = self._find_newcps(cps_low)  
        self._new_cps_up = self._find_newcps(cps_up)
        return self._new_cps_low, self._new_cps_up
        
    def _find_newcps(self,cps):
        """
        Finds the new curve control points dependent on the amount of 
        iterations desired.
        """
        self._it = it        
        
        # convert points list to normal list
        A = list(zeros(len(cps)))
        for i in xrange(0,len(A)):
            A[i] = cps[i].x,cps[i].y,cps[i].z            
        
        n = 0
        
        # if no iteration required simply return the control points
        if it <= 0:
                n = it + 1
                return cps            
        
        # start loop to find new control points
        while n < it:
            
            n += 1
#            print "n = %s -------------- " %(n)           
            
            # find intermediate points 
            xmid = list(zeros(len(A)-3))
            ymid = list(zeros(len(A)-3))
            zmid = list(zeros(len(A)-3))    
                
            for i in xrange(0,2**(n-1)-1):
                xmid[i] = A[i+1][0]*0.5 + A[i+2][0]*0.5
                ymid[i] = A[i+1][1]*0.5 + A[i+2][1]*0.5
                zmid[i] = A[i+1][2]*0.5 + A[i+2][2]*0.5
#            print 'xmid = %s' %xmid
#            print 'ymid = %s' %ymid
#            print 'A = %s' %A
            
            # contruct new list with points and intermediate coordinates 
            xtot = list(zeros(len(A)+len(xmid)))
            ytot = list(zeros(len(A)+len(xmid)))
            ztot = list(zeros(len(A)+len(xmid)))
            
            for i in range(0,2):
                    xtot[i] = A[i][0]
                    ytot[i] = A[i][1]
                    ztot[i] = A[i][2]
                    
            for i in range(2,len(xtot)-2):
                    if i % 2 == 0:
                        xtot[i] = xmid[i/2-1]
                        ytot[i] = ymid[i/2-1]
                        ztot[i] = zmid[i/2-1]
                    if i % 2 != 0:
                        xtot[i] = A[int(i/2)+1][0]
                        ytot[i] = A[int(i/2)+1][1]
                        ztot[i] = A[int(i/2)+1][2] 
            
            if n == 1:                    
                for i in range(len(A)-1,len(A)):
                    xtot[i+len(xmid)] = A[i][0]
                    ytot[i+len(xmid)] = A[i][1]
                    ztot[i+len(xmid)] = A[i][2]
                                      
            else:
                for i in range(len(A)-2,len(A)):
                    xtot[i+len(xmid)] = A[i][0]
                    ytot[i+len(xmid)] = A[i][1]
                    ztot[i+len(xmid)] = A[i][2]
                
#            print 'xtot = %s' %xtot   
#            print 'ytot = %s' %ytot  
            
            # construct the final cooridinates
            xfin = list(zeros(len(xtot)+1))
            yfin = list(zeros(len(xtot)+1))
            zfin = list(zeros(len(xtot)+1))
            
            xfin[0] = xtot[0]
            yfin[0] = ytot[0]
            zfin[0] = ztot[0]
                    
            for i in range(len(xtot)-1,len(xtot)):
                xfin[i+1] = xtot[i]
                yfin[i+1] = ytot[i]
                zfin[i+1] = ztot[i]
                
            for i in range(1, len(xtot)):
                xfin[i] = xtot[i-1]*0.5 + xtot[i]*0.5
                yfin[i] = ytot[i-1]*0.5 + ytot[i]*0.5
                zfin[i] = ztot[i-1]*0.5 + ztot[i]*0.5
                      
#            print 'xfin = %s' %xfin
#            print 'yfin = %s' %yfin
             
            # construct new list of control points 
            C = list(zeros(2+2**n))
            for i in xrange(0,len(C)):
                C[i] = xfin[i],yfin[i],zfin[i]
                
#            print 'C = %s' %C
            A = C[:]
            self._A = A     
            
            # Make the points list
            Apoints = A[:]
            for i in range(0,len(A)):
                Apoints[i] = Point(A[i][0],A[i][1],A[i][2])
            self._Apoints = Apoints
#            print Apoints
            
        return self._Apoints
    
        
    
    
if __name__=='__main__':
    
    close('all')

    test = MeridionalChannel(height,length,R1,R2,outlet_angle,it)
    testspline1 = Bspline(test._cps_low)
    testspline2 = Bspline(test._cps_up)
    
    x1 = zeros(len(test._cps_low))
    y1 = list(x1)
    x2 = list(x1)
    y2 = list(x1)
    
    for i in range(0,len(test._cps_low)):
        x1[i] = test._cps_low[i].x
        y1[i] = test._cps_low[i].y
        
    for i in range(0,len(test._cps_up)):
        x2[i] = test._cps_up[i].x
        y2[i] = test._cps_up[i].y
    
    fig1 = plt.figure()
    title('Meridional Channel')
    axis([-0.6,length*1.1,-0.5,height*1.1])  
    testspline1.plot()
    testspline2.plot()
    scatter(x1,y1,label = 'Control Points')
    scatter(x2,y2)
    plt.plot([test._cps_low[0].x,test._cps_up[0].x],[test._cps_low[0].y,test._cps_up[0].y],'g',label = 'Inlet')
    plt.plot([test._cps_low[-1].x,test._cps_up[-1].x],[test._cps_low[-1].y,test._cps_up[-1].y],'r',label = 'Outlet')
#    plt.plot([test._cps_low[-1].x,0.0],[test._cps_low[-1].y,-length*tan(outlet_angle*pi/180.0)],'y-')
    xlabel('length x')
    ylabel('height y')
    legend()

        
        
        
        
        
        
        
        
        
        