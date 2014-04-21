"""
Created on Apr 18, 2014

@author: salvatorevitale
"""

from point import Point, Norm, Angle
from math import sin, cos
from matplotlib.pylab import *

class Line(object):
    """
    create a parametric line  
    """


    def __init__(self, P, T):
        """
        to costruct a line it needs a point and a direction
        P -> Point
        T -> norm
        
        """
        self.P = P
        self.T = T
        
    def __call__(self, u):
        P, T = self.P, self.T
        return Point(P.x + T.x*u, P.y + T.y*u, P.z + T.z*u)
    
    def intersect(self, other):
        P1,T1,P2,T2 = self.P, self.T, other.P, other.T
        u = (T2.y*(P2.x - P1.x) - T2.x*(P2.y - P2.y))/(T1.x*T2.y - T1.y*T2.x)
        return self(u)
           
    def plot(self):
        P1 = self.P
        P2 = self.__call__(1.0)
        plot([P1.x, P2.x],[P1.y, P2.y])
        
    
    
    
    
    
if __name__=='__main__':
    
    ang1 = Angle(30.0)
    ang2 = Angle(-60.0)
    P1 = Point()
    P2 = Point(1.0,0.0)
    T1 = Norm(cos(ang1()), sin(ang1()))
    T2 = Norm(cos(ang2()), sin(ang2()))
    l1 = Line(P1,T1)
    l2 = Line(P2,T2)
#     l1.plot()
#     l2.plot()
    T1.plot2D()
    print(l1.intersect(l2))
    show()
    
  