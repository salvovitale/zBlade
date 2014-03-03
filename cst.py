import sys
utility_folder = '../utility/'
sys.path.insert(0, utility_folder)
import csv
from math import *
from matplotlib.pylab import *
from sys import argv
import util as ut
import preproc_prof as pp
import scipy.optimize



#class Profile:
      



class Cst:
    
    def __init__(self, a):
        self._N = a[0:2]
        self._A = a[2:len(a)]
    
    def __call__(self, x):
        n = len(self._A)
        S = 0
        for i in xrange(n):
            S += self._A[i]*self._bin_coef(i,n-1)*(1-x)**(n-1-i)*x**i
        return self._class_function(x)*S

    def _class_function(self, x):
        return x**self._N[0]*(1-x)**self._N[1]

    def _bin_coef(self, i, n):
        return factorial(n)/(factorial(i)*factorial(n-i))

        
class BF_cst:
    
    def __init__(self, xy_p, c_err = 1.0, n_ord = 2.0):
        self._xy_p = xy_p
        self._c_err = c_err
        self._n_ord = n_ord

    def __call__(self, A):
        self._cst = Cst(A)
        return np.linalg.norm(self._err(), ord = self._n_ord)
        
    def _err(self):
        return np.where(self._xy_p[1,:] != 0.0, abs(self._xy_p[1,:] - self._cst(self._xy_p[0,:]))/abs(self._xy_p[1,:])**self._c_err, abs(self._xy_p[1,:] - self._cst(self._xy_p[0,:])))
        


    @staticmethod
    def init_cfc( cfc_ig1 = 0.4, cfc_ig2 = 0.4):
        N = np.ones(2)
        N[0] *= cfc_ig1
        N[1] *= cfc_ig2
        return N

    @staticmethod
    def init_A(poly_gr = 5, A_ig = 0.8):
        A = np.ones(poly_gr+1)
        A *= A_ig
        return A

    @staticmethod    
    def bound_cfc( cfc_lb1 = 0.1, cfc_ub1 = 1.5, cfc_lb2 = 0.1 , cfc_ub2 = 1.4):
        N = np.ones((2,2), dtype = float)
        N[0,0] *= cfc_lb1
        N[1,0] *= cfc_lb2
        N[0,1] *= cfc_ub1
        N[1,1] *= cfc_ub2
        return N
        
    @staticmethod
    def bound_A(poly_gr = 4, A_lb = -1.5, A_ub = +1.5):
        A  = np.ones((poly_gr + 1, 2), dtype = float)
        A[:,0] *= A_lb
        A[:,1] *= A_ub
        return A


        




if __name__=='__main__':
    script, input_file = argv
    
    import ConfigParser
    config = ConfigParser.ConfigParser()
    config.read(input_file)
    filename = config.get("general", "filename")                 
    xy_ss, xy_ps = pp.import_prof(filename)

    poly_gr_ss = int(config.get("ss", "poly_gr_ss"))
    cfc_lb1_ss = float(config.get("ss", "cfc_lb1"))
    cfc_ub1_ss = float(config.get("ss", "cfc_ub1"))
    cfc_lb2_ss = float(config.get("ss", "cfc_lb2"))
    cfc_ub2_ss = float(config.get("ss", "cfc_ub2"))
    cfc_ig1_ss = float(config.get("ss", "cfc_ig1"))
    cfc_ig2_ss = float(config.get("ss", "cfc_ig2"))
    A_lb_ss = float(config.get("ss", "A_lb"))
    A_ub_ss = float(config.get("ss", "A_ub"))    
    A_ig_ss = float(config.get("ss", "A_ig"))
    c_ss    = float(config.get("ss", "yp_c"))

    poly_gr_ps = int(config.get("ps", "poly_gr_ps"))
    cfc_lb1_ps = float(config.get("ps", "cfc_lb1"))
    cfc_ub1_ps = float(config.get("ps", "cfc_ub1"))
    cfc_lb2_ps = float(config.get("ps", "cfc_lb2"))
    cfc_ub2_ps = float(config.get("ps", "cfc_ub2"))
    cfc_ig1_ps = float(config.get("ps", "cfc_ig1"))
    cfc_ig2_ps = float(config.get("ps", "cfc_ig2"))
    A_lb_ps = float(config.get("ps", "A_lb"))
    A_ub_ps = float(config.get("ps", "A_ub"))    
    A_ig_ps = float(config.get("ps", "A_ig"))
    c_ps    = float(config.get("ps", "yp_c"))
    
    N_ss  = BF_cst.init_cfc( cfc_ig1 = cfc_ig1_ss, cfc_ig2 = cfc_ig2_ss)
    N_ps  = BF_cst.init_cfc( cfc_ig1 = cfc_ig1_ps, cfc_ig2 = cfc_ig2_ps)
    A_ss  = BF_cst.init_A(poly_gr = poly_gr_ss, A_ig = A_ig_ss)
    A_ps  = BF_cst.init_A(poly_gr = poly_gr_ps, A_ig = A_ig_ps)


    Ass_o = np.concatenate((N_ss, A_ss))
    Aps_o = np.concatenate((N_ps, A_ps))
    
    Nss_b = BF_cst.bound_cfc( cfc_lb1 = cfc_lb1_ss, cfc_ub1 = cfc_ub1_ss, cfc_lb2 = cfc_lb2_ss, cfc_ub2 = cfc_ub2_ss)
    Nps_b = BF_cst.bound_cfc( cfc_lb1 = cfc_lb1_ps, cfc_ub1 = cfc_ub1_ps, cfc_lb2 = cfc_lb2_ps, cfc_ub2 = cfc_ub2_ps)
    Ass_b = BF_cst.bound_A(poly_gr = poly_gr_ss, A_lb = A_lb_ss, A_ub = A_ub_ss)
    Aps_b = BF_cst.bound_A(poly_gr = poly_gr_ps, A_lb = A_lb_ps, A_ub = A_ub_ps)

    Ass_b = np.vstack((Nss_b, Ass_b))
    Aps_b = np.vstack((Nps_b, Aps_b))

    of_s = BF_cst(xy_ss, c_err = c_ss)
    of_p = BF_cst(xy_ps, c_err = c_ps)
    

    
#    Ass, nfeval_ss, rc_ss = scipy.optimize.fmin_tnc(of_ss, Ass_o, approx_grad=True, bounds = Ass_b, maxfun = 2000)
#    Aps, nfeval_ps, rc_ps = scipy.optimize.fmin_tnc(of_ps, Aps_o, approx_grad=True, bounds = Aps_b, maxfun = 2000)


    Ass = scipy.optimize.fmin_slsqp(of_s, Ass_o, bounds = Ass_b, iter = 1000)
    Aps = scipy.optimize.fmin_slsqp(of_p, Aps_o, bounds = Aps_b, iter = 1000)

#    print scipy.optimize.anneal(of_ss, Ass_o, maxeval = 2000)
#    print Ass, Aps
#    Ass = np.array([0.5505, 0.3403, 0.6500, 0.5169, 0.1713, 0.2214, 0.0375, 0.0397, 0.0108, 0.0173])
#    Aps = np.array([0.4612, 0.2358, -0.4779, -1.0122, 0.1103, -0.0936, -0.0266, 0.0121, 0.0008, - 0.0101])

    print Aps, Ass

#    cst_ps = cst(xy_ps[0,:], Aps)
#    cst_ss = cst(xy_ss[0,:], Ass)
    cst_ps = Cst(Aps)
    cst_ss = Cst(Ass)

   # print of_ss(Ass), of_ps(Aps)
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex= True, sharey = True)
    
    ax1.plot(xy_ps[0,:], xy_ps[1,:] , 'r-', xy_ss[0,:],xy_ss[1,:], 'r-')
    ax1.plot(xy_ss[0,:], cst_ss(xy_ss[0,:]), 'b-', xy_ps[0,:], cst_ps(xy_ps[0,:]), 'b-')
#    ax1.plot(xy_ss[0,:], cst_ss(xy_ss[0,:]), 'b-', xy_ps[0,:], cst_ps(xy_ps[0,:]), 'b-')
    ax1.axis('equal')
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
   # ax1.set_legend(['real profile'], ['paramatrized profile'])

    p21, p22, p23, p24, = ax2.plot(xy_ps[0,:], xy_ps[1,:], 'r-', xy_ps[0,:], xy_ps[1,:], 'bo', xy_ss[0,:],xy_ss[1,:],'r-', xy_ss[0,:], xy_ss[1,:], 'bo')
    p25,  = ax2.plot(xy_ss[0,:], cst_ss(xy_ss[0,:]), 'b-')
    p26,  = ax2.plot(xy_ps[0,:], cst_ps(xy_ps[0,:]), 'b-')
    legend([p21, p24, p26], ["real", "real points", "parametrized"], bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
   # plt.plot(xy_ps[0,:], xy_ps[1,:], 'bo')
   # plt.plot(xy_ss[0,:], xy_ss[1,:])
   # plt.plot(xy_ss[0,:], xy_ss[1,:], 'bo')
    
    ax3.plot(xy_ps[0,:], xy_ps[1,:], 'r-')
    ax3.plot(xy_ss[0,:], xy_ss[1,:], 'r-')
    
    ax4.plot(xy_ss[0,:], cst_ss(xy_ss[0,:]), 'b-')
    ax4.plot(xy_ps[0,:], cst_ps(xy_ps[0,:]), 'b-')
    plt.show()
    raw_input('\>')


#plot(x,y)
#show()


