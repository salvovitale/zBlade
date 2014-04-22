import sys
 
utility_folder = '../utility/'
sys.path.insert(0, utility_folder)
from matplotlib.pylab import *
from math import factorial
import util as ut
import preproc_prof as pp
import scipy.optimize
from point import Point
from plot import Plot
from bezier import Bezier






    
    
    
class Basis_bspline:
        
    def __init__(self, U, p = 1):
        """
        construct Bspline basis function object

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
        self._i = self._find_span(u)
        i = self._i
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
             
        return mid
        
    def get_span_index(self):
        return self._i    




class Bspline:
    """
    This class implement a rational Bspline with a nonperiodic uniform
    knot vectors (definition at pag 66-67 of Nurbs book). 
    """
    def __init__(self, P, p = 2, a = 0, b = 1, kth = 0):
        """
        Constructor 
        input:
        p == degree of the curve 
        P == list of control point   
        a == lower bound of the interval for the knot vector domain
        b == upper bound of the interval for the knot vector domain
        
        
        identities:
        n+1 = len(P) number of cont. Point namely number of basis function
        m+1 = len(U) length of the knot vector
        m = n + p + 1
        """
        self._P = P
        self._p = p
        self._a = a
        self._b = b
        self._n = len(P) - 1
        self._m = self._n + p + 1
        self._kth = kth
        self._U = self._knotsvector()
        self._X, self._Y, self._Z, self._W = self._sep()
        
    def __call__(self, u):
        """
        this method evaluats the bspline rational function at the value of the parameter u 
        see pag 82 of Nurbs book
        
        """
        p, U, kth, X, Y, Z, W =  self._p, self._U, self._kth, self._X, self._Y, self._Z, self._W  
        N = DerBasisBspline(U = U, p = p)
        N(u)
        span = N.get_span_index()
        x = sum(N(u)[kth][:]*X[span - p : span+1]*W[span - p : span+1])
        y = sum(N(u)[kth][:]*Y[span - p : span+1]*W[span - p : span+1])
        z = sum(N(u)[kth][:]*Z[span - p : span+1]*W[span - p : span+1])
        w = sum(N(u)[kth][:]*W[span - p : span+1])
        return x/w, y/w, z/w
    
    def _knotsvector(self):
        
        """
        construct a uniform nonperiodic knots vectors
        {a_0, a_1,..., a_p , u_p+1,..., u_m-p-1, b_m-p, ... ,b_m-1, b_m}
        
        d = (b-a)/(m-2p)
        u_p+1 = a + d
        u_p+2 = a + 2*d
        ...
        u_m-p-1 = a + j*d 
        """
        p, a, b = self._p, self._a, self._b
        U1 = np.ones(p+1)
        U1 *= a
        U2 = self._internal_knots() 
        U3 = np.ones(p+1)
        U3 *= b
        return np.concatenate((U1, np.concatenate((U2, U3))))
                      

    def _internal_knots(self):
        """
        calculate the values of the uniform internal knots
          
        """
        p, m, a, b = self._p, self._m, self._a, self._b
        d = (b-a)/float(m-2*p)
        int_knot = np.ones((m-2*p-1))
        for j in xrange(m-2*p-1):
            int_knot[j] = a + (j+1)*d
        return int_knot
    
    def _sep(self):
        n, P = self._n, self._P
        x = np.ones(n+1)
        y = np.ones(n+1)
        z = np.ones(n+1)
        w = np.ones(n+1)
        for i in xrange(n+1):
            x[i], y[i], z[i], w[i] = P[i].split()
             
        return x, y, z, w
    
    def get_int_U(self):
        p , m = self._p, self._m
        return self._U[p: m-p+1]
    
    def get_U(self):
        return self._U
    
    def get_n(self):
        return self._n
    
    def get_p(self):
        return self._p
    
    def get_m(self):
        return self._m
    
    def get_cp(self):
        return self._P
    
    def get_ncp(self):
        return self._n + 1
        
    def get_x(self):
        return self._X
    
    
    def get_y(self):
        return self._Y
    
    def get_z(self):
        return self._Z
    
    
    def get_w(self):
        return self._W
    
    def plot(self):
        a = Plot(self)
        a()






class DerBasisBspline:
        
    def __init__(self, U, p = 2):
        """
        construct Bspline basis function object

        p == degree of the curve 
        U == nonperiodic Knot vector  
        nd == max derivative order
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
#         if nd <= p:
#             self._nd = nd
#         else:
#             print "you cannot calcualte derivative for k major than p"
        self._nd = p
                
    def __call__(self, u):
        """
        this method compute the non vanishing basis function and their derivatives
        Nurb book pag 72-73 using Eq. 2.9 pag 61 
        """
        self._i = self._find_span(u)
        i = self._i
        ders = self._compute_basisFunction(u)
        return ders
                
    def _compute_basisFunction(self, u):
        p, U, i, nd  = self._p, self._U, self._i, self._nd
        N = np.ones((p+1, p+1), dtype = float)
        left = np.ones(p+1)
        right = np.ones(p+1)
        for j in xrange(1, p+1, 1):
            left[j] = u - U[i+1-j]
            right[j] = U[i+j]- u
            saved = 0.0
            for r in xrange(j):
                #lower triangle
                N[j][r] = right[r+1] + left[j-r]
                temp = N[r][j-1]/N[j][r]
                #upper triangle
                N[r][j] = saved + right[r+1]*temp
                saved = left[j-r]*temp
            N[j][j] = saved    
        ders = np.ones((nd+1, p+1), dtype = float) 
        a = np.ones((2, p+1), dtype = float)
        #load the basis functions
        for j in xrange(p+1):
            ders[0][j] = N[j][p]
        # This section computes the derivatives using Eq. 2.9 
        for r in xrange(p+1):
            s1 = 0
            s2 = 1
            a[0][0] = 1.0
            #loop to compute kth derivatives
            for k in xrange(1, nd+1, 1): 
                d = 0.0
                rk = r - k
                pk = p - k
                if r >= k:
                    a[s2][0] = a[s1][0]/N[pk+1][rk]
                    d = a[s2][0]*N[rk][pk]
                if rk >= -1:
                    j1 = 1
                else:
                    j1 = -rk
                if r-1 <= pk:
                    j2 = k-1
                else:
                    j2 = p-r
                for j in xrange(j1, j2+1,1):
                    a[s2][j] = (a[s1][j] - a[s1][j-1])/N[pk+1][rk+j]
                    d += a[s2][j]*N[rk+j][pk]
                if r <= pk:
                    a[s2][k] = -a[s1][k-1]/N[pk+1][r]
                    d += a[s2][k]*N[r][pk] 
                ders[k][r] = d
                j = s1
                s1 = s2
                s2 = j
            #multiply through by the correct factor
        r = p
        for k in xrange(1, nd+1, 1):
            for j in xrange(p+1):
                ders[k][j] *= r
                r *= (p-k)
        return ders         
                    
     
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
             
        return mid
        
    def get_span_index(self):
        return self._i    
    
    
    
    

            
            
if __name__=='__main__':
#     U = np.array([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0])
#     p = 2         
#     u = 2.5
#     prova = DerBasisBspline(U)
#     print prova(2.5)[0][2]
#     print prova(2.5)[1][2]
#     print prova(2.5)[2][2]   
#     prova = Basis_bspline(U,p)




    P = [Point(0.0, 0.0), Point(0.5, 1.0), Point(1.0, 0.0), Point(1.5, 2.0), Point(3.5, -1.0), Point(4.5, 0.0) ]
    p = 5
    prova1 = Bspline(P, p, kth = 3)
    prova2 = Bspline(P, p, kth = 2)
    prova1.plot()
    prova2.plot()
    axis('equal') 
    show()  
# #     P = [Point(0.0, 1.0), Point(2.79761904762, 1.59530618983), Point(5.59523809524, 2.2824786079), Point(8.39285714286, 2.50866465597), Point(11.1904761905, 2.86509988565)]
# #     print prova(4.0)
# #     print sum(prova(4.0))
#     prova2 = Bspline(P, p) 
#     print prova2.__class__.__name__   
#     prova2.plot()
#     print prova2.get_int_U() 
#     prova3 = Bezier(P)
#     prova3.plot()
#     axis('equal') 
#     show()      
    
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