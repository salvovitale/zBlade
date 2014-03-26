import sys
utility_folder = '../utility/'
sys.path.insert(0, utility_folder)
from matplotlib.pylab import *
from math import factorial
import util as ut
import preproc_prof as pp
import scipy.optimize





class Point:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0 , w = 1.0 ):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        
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
        return "Point(%s, %s)" % (self.x, self.y)      
    
    def split(self):
        return self.x, self.y, self.z, self.w
    
    
    
class Basis_bspline:
        
    def __init__(self, U, p = 0):
        """
        construct Bspline basis function object
        uses Cox-DeBoor

        p == degree of the curve 
        U == nonperiodic Knot vector  
        
        identities:
        m = len(U) - 1
        n = m - p - 1
        
        m+1 == length of the knot vector
        n+1 == number basis function
        """
        self._U = U
        self._p = p
        self._m = len(U) - 1
        self._n = self._m - p - 1
            
    def __call__(self, u):
        """
        this method compute the non vanishing basis function
        Nurb book pag 70
        """
        p, U = self._p, self._U
        self._find_span(u)
        i = self.get_span_index()
        N = np.ones(p+1)
        left = np.ones(p+1)
        right = np.ones(p+1)
        for j in xrange(1, p+1, 1):
            left[j] = u - U[i+1-j]
            right[j] = U[i+j]- u
            saved = 0.0
            for r in xrange(j):
                temp = N[r]/(right[r+1] + left[j-r])
                N[r] = saved + right[r+1]*temp
                saved = left[j-r]*temp
            N[j] = saved
        return N
                
        
        
        return u 
     
    def _find_span(self, u):
        """
        this method determine the knot span index
        Nurbs book pag 68
        
        first it handles a special case u == um (last element of the knot vector)
        then it uses a binary search
        
        """
        n, p, U = self._n, self._p, self._U
        
        if u == U[n+1]:
            return n
        low = p
        high = n+1
        mid = (low + high)/2
        while u < U[mid] or u >= U[mid +1]:
            if u < U[mid]:
                high = mid
            else:
                low = mid
            mid = (low+high)/2    
        self._i = mid
        
    def get_span_index(self):
        return self._i    



if __name__=='__main__':
    U = np.array([0, 0, 0, 1, 2, 3, 4, 4, 5, 5, 5])
    p = 3            
    prova = Basis_bspline(U,p)
    print prova(3.0)
    print sum(prova(3.0))        
# class Bspline(object):
# 
#     def __init__(self, P, u, k = None):
#         """
#         construct Bspline object
#         uses Cox-DeBoor
# 
#         P == vector of two-dimensional control points
#         t == vector of non-decreasing real numbers
#         k == degree of curve
# 
#         identities:
#         P = (P[0], ... P[n]); n = len(P) - 1
#         t = (t[0], ... t[m]); m = len(t) - 1
#         k = m - n - 1
#         m = n + k + 1
#         n = m - k - 1
#         """
#         m, n = len(t) - 1, len(P) - 1
#         if not k: k = m - n - l
#         else: assert m == n + k + 1
#         self.k, self.t = k, t
#         self.X, self.Y = _unzip(P) # points in X, Y components
#         self._deboor() # evaluate
# 
#     def __call__(self, t_):
#         """
#         S(t) = sum(b[i][k](t) * P[i] for i in xrange(0, n))
#         domain: t in [t[k - 1], t[n + 1]]
# 
#         returns point on Bspline at t_
#         """
#         k, t = self.k, self.t
#         m = len(t) - 1
#         n = m - k - 1
#         assert t[k - 1] <= t_ <= t[n + 1] # t in [t[k - 1], t[n + 1]]
#         X, Y, b = self.X, self.Y, self.b
#         x, y, _n = 0, 0, xrange(n + 1) # initial return values, iterator over P
#         for i in _n:
#             b_i = b[i][k](t_)
#             x += X[i] * b_i
#             y += Y[i] * b_i
#         return x, y
# 
#     def _deboor(self):
#         # de Boor recursive algorithm
#         # S(t) = sum(b[i][k](t) * P[i] for i in xrange(0, n))
#         #
#         # b[i][k] = {
#         #     if k == 0:
#         #         t[i] <= t_ < t[i+1]
#         #     else:
#         #         a[i][k](t)*b[i][k-1](t)+(1-a[i+1][k](t))*b[i+1][k-1](t)
#         # }
#         #
#         # a[i][k] = {
#         #     if t[i] == t[i+k]:
#         #         0
#         #     else:
#         #         (t_-t[i])/(t[i+k]-t[i])
#         # }
#         #
#         # NOTE: for b[i][k](t), must iterate to t[:-1];
#         # the number of [i, i + 1) spans in t
#         k, t = self.k, self.t
#         m = len(t) - 1 # iterate to t[:-1]
#         a, b, _k_, _m_ = [], [], xrange(k + 1), xrange(m)
#         for i in _m_:
#             a.append([]); b.append([]) # a[i]; b[i]
#             for k in _k_:
#                 a[i].append(None) # a[i][k]
#                 # if k == 0: b[i][k](t) is a step function in [t[i], t[i + 1])
#                 if k == 0: b[i].append(lambda t_, i=i: t[i] <= t_ < t[i + 1])
#                 # if m < i + k: b[i][k](t) undefined
#                 elif m < i + k: b[i].append(lambda t_: False)
#                 # else: calculate b[i][k](t)
#                 else:
#                     # if t[i] == t[i + k]: a[i][k] undefined
#                     if t[i] == t[i + k]: a[i][k] = lambda t_: False
#                     # else: calculate a[i][k](t)
#                     else:
#                         # a[i][k](t) = (t_ - t[i]) / (t[i + k] - t[i])
#                         a[i][k] = lambda t_, i=i, k=k: ((t_ - t[i]) /
#                                                         (t[i + k] - t[i]))
#                     # b[i][k](t) = a[i][k](t) * b[i][k - 1](t) +
#                     #              (1 - a[i + 1][k](t)) * b[i + 1][k - 1](t)
#                     b[i].append(lambda t_, i=i, k=k:
#                                 a[i][k](t_) * b[i][k - 1](t_) +
#                                 (1 - a[i + 1][k](t_)) * b[i + 1][k - 1](t_))
#         self.b = b
# 
#     def insert(self, t_):
#         """
#         Q[i] = (1 - a[i][k]) * P[i] + a[i][k] * P[i]
#         domain: t in (t[0], t[m])
# 
#         insert new control point at t_
#         """
#         t = self.t
#         assert t[0] < t_ < t[-1] # t_ in (t[0], t[m])
#         X, Y, k = self.X, self.Y, self.k
#         m = len(t) - 1
#         _t_ = xrange(m + 1)
#         # find the span containing t_
#         for i in _t_:
#             if t[i] <= t_ < t[i + 1]: break
#         assert not i < k + 1 and not i > m - k + 1 # i not in clamp
#         Q_x, Q_y = [], [] # new control points
#         # iterate over replaced control points
#         # set new control points
#         for j in xrange(i - k + 1, i + 1):
#             a_j = (t_ - t[j]) / (t[j + k] - t[j])
#             Q_x.append((1 - a_j) * X[j - 1] + a_j * X[j])
#             Q_y.append((1 - a_j) * Y[j - 1] + a_j * Y[j])
#         Q_x, Q_y = tuple(Q_x), tuple(Q_y)
#         self.t = t[:i + 1] + [t_] + t[i + 1:]
#         self.X = X[:i - k + 1] + Q_x + X[i:]
#         self.Y = Y[:i - k + 1] + Q_y + Y[i:]
#         self._deboor() # re-evaluate