from matplotlib.pylab import *
from point import Point
from bspline import Bspline





class NeedleDef(object):
    """
    we define the needle to be inserted into the  divergent part of the nozzle
    
    
    """


    def __init__(self, tip_pos = 2.0, del_D = 2.0, alpha = 30.0, r_curv = 0.2, r_curv_c = 0.2, length = 5.0):
        """
        Constructor
        
        l-> lenght
        A-> inlet Area
        m-> inlet derivative
        ncp-> number of Control Point
        type -> stretching function
        """
        alpha_rad = 3.1415*alpha/180.0
        self._tip_pos = tip_pos
        self._del_D = del_D
        self._alpha = alpha_rad
        self._r_curv = r_curv
        self._r_curv_c = r_curv_c
        self._length = length
        self._cpt, self._cpc, self._cps = self._define_cp()
        
   
   
    def _define_cp(self):
        tip_pos, del_D, alpha, r_curv, r_curv_c, length = self._tip_pos, self._del_D, self._alpha, self._r_curv, self._r_curv_c, self._length
        P0 = Point(tip_pos, 0.0)
        y1 = r_curv*del_D
        P1 = Point(tip_pos, y1)
        del_x = del_D/tan(alpha)
        P2f = Point(tip_pos + del_x , del_D)
        del_xc = del_x*r_curv_c 
        del_yc = del_xc*(del_D - y1)/del_x
        
        P2 = Point(tip_pos + del_x - del_xc , del_D - del_yc)
        P3 = Point(tip_pos + del_x + del_xc , del_D)
        P4 = Point(tip_pos + length , del_D)
        return [P0, P1, P2],[P2, P2f, P3], [P3, P4]
         
        
    def getCPt(self):
        return self._cpt
    
    def getCPc(self):
        return self._cpc
    
    def getCPs(self):
        return self._cps
     


class Tip(object):
    """
    we define the needle to be inserted into the  divergent part of the nozzle
    
    
    """


    def __init__(self, needleDef):
        """
        Constructor
        
        l-> lenght
        A-> inlet Area
        m-> inlet derivative
        ncp-> number of Control Point
        type -> stretching function
        """
        self._needleDef = needleDef
        self._P = needleDef.getCPt()
        self._tip = Bspline(self._P, p = 2)
                
    def __call__(self, u):
        return  self._tip(u)          
        
        
    def plot(self):
        return self._tip.plot()  
    

class Corner(object):
    """
    we define the needle to be inserted into the  divergent part of the nozzle
    
    
    """


    def __init__(self, needleDef):
        """
        Constructor
        
        l-> lenght
        A-> inlet Area
        m-> inlet derivative
        ncp-> number of Control Point
        type -> stretching function
        """
        self._needleDef = needleDef
        self._P = needleDef.getCPc()
        self._corner = Bspline(self._P, p = 2)
                
    def __call__(self, u):
        return  self._corner(u)          
        
        
    def plot(self):
        return self._corner.plot()        
    
    
class Shaft(object):
    """
    we define the needle to be inserted into the  divergent part of the nozzle
    
    
    """


    def __init__(self, needleDef):
        """
        Constructor
        
        l-> lenght
        A-> inlet Area
        m-> inlet derivative
        ncp-> number of Control Point
        type -> stretching function
        """
        self._needleDef = needleDef
        self._P = needleDef.getCPs()
        self._shaft =  Bspline(self._P, p = 1)
        
    def __call__(self, u):
        return  self._shaft(u)          
        
        
    def plot(self):
        return self._shaft.plot()    
        
        
        