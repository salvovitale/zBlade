""" @package zBlade
In this module it is implemented the class Point which represent 
the cartesian point    
"""

from math import sqrt, pi, pow
from matplotlib.pylab import plot


class Point:
    """ class Point
    
    It is implemented as a 4D point in order to generalize
    the calculation of rational curves. 
    See Bspline and Bezier to understand what this means   
        
    """
    def __init__(self, x = 0.0, y = 0.0, z = 0.0 , w = 1.0 ):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        
    def __call__(self):
        return self.x, self.y, self.z     
        
    def distance(self, other):
        return sqrt((other.x-self.x)*(other.x-self.x)+(other.y-self.y)*(other.y-self.y)+(other.z-self.z)*(other.z-self.z))

    def length(self):
        return self.distance(Point(0.0, 0.0, 0.0))

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y, self.z-other.z)

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y, self.z + other.z)

    def __mul__(self, c):
        return Point(c*self.x, c*self.y, c*self.z)

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other):
        return not (self == other)
    
    def towards(self, target, t):
        return Point((1.0-t)*self.x+t*target.x, (1.0-t)*self.y+t*target.y, (1.0-t)*self.z+t*target.z)
    
    def halfway(self, target):
        return Point((self.x+target.x).div2(), (self.y+target.y).div2(), (self.z+target.z).div2())

    def compare_lex(self, other):
        if self.x < other.x:
            return -1
        if self.x > other.x:
            return 1
        if self.y < other.y:
            return -1
        if self.y > other.y:
            return 1
        return 0

    def less_lex(self, p):
        return self.compare_lex(p) < 0

    def less_eq_lex(self, p):
        return self.compare_lex(p) <= 0
    
    def __repr__(self):
        return "Point(%s, %s, %s, %s)" % (self.x, self.y, self.z, self.w)      
    
    def split(self):
        return self.x, self.y, self.z, self.w



class Vector:
    
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z
        

    def __repr__(self):
        return "Vector(%s, %s, %s)" % (self.x, self.y, self.z) 
    
    def norm(self):
        return normVector(Vector(self.x, self.y, self.z))



# class Norm:
#     
#     def __init__(self, v):
#         self.x = x/sqrt(x**2+y**2+z**2)
#         self.y = y/sqrt(x**2+y**2+z**2)
#         self.z = z/sqrt(x**2+y**2+z**2)
#     
#     def __call__(self):
#         return self.x, self.y, self.z#, [0.0, self.z]
#     
#     def plot2D(self):
#         plot([0.0, self.x], [0.0, self.y])
        
       
            
        
class Angle:
    """
    to convert angle from degree to radiant
    """
    def __init__(self, angle):
        self._angle = angle
    def __call__(self):
        return pi*self._angle/180.00 
    def getAngle(self):
        return self._angle          
    

def normVector(v):
    x = v.x/sqrt(v.x**2+v.y**2+v.z**2)
    y = v.y/sqrt(v.x**2+v.y**2+v.z**2)
    z = v.z/sqrt(v.x**2+v.y**2+v.z**2)
    return Vector(x, y, z)
            
def normal2D(v):
    return Vector(-v.y, v.x)