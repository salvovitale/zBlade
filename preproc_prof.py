import sys
utility_folder = '../utility/'
sys.path.insert(0, utility_folder)
import csv
from math import *
from matplotlib.pylab import *
from sys import argv
import util as ut



def distance_2p(xy1,xy2):
    return  np.sqrt((xy1[0] - xy2[0])**2 + (xy1[1] - xy2[1])**2)

def find_le(xy_p,xy_te):
    dist = (xy_p[0,:] - xy_te[0])**2 + (xy_p[1,:] - xy_te[1])**2
    return dist.argmax()

#def shift_0

def read_xy_p (filename):
    data = ut.read_data(filename)
    return np.array([[row[0] for row in data], [row[1] for row in data]])


def shift_o(xy_p, xy_le):

    xy_p[0,:] -= xy_le[0]
    xy_p[1,:] -= xy_le[1]
    return xy_p


def atang(xy):
    den = xy[0]**(-1)
    alpha = np.arctan(xy[1]*den)*180.0/pi
    import operator
    condition1 = operator.and_(xy[0] < 0.00, xy[1] < 0.0)
    condition2 = operator.and_(xy[0] > 0.0, xy[1] > 0.0)
    condition3 = operator.and_(xy[0] > 0.0, xy[1] < 0.0)
    condition4 = operator.and_(xy[0] < 0.0, xy[1] > 0.0)
    condition23 =operator.or_(condition2, condition3)
    r = np.where(condition4, alpha + 180.0,
        np.where(condition23, alpha, 
        np.where(condition1, alpha - 180.0, 0.0)))

    return r


def rotate(xy_p, ang):
    
    
    r = sqrt(xy_p[0,:]**2 + xy_p[1,:]**2)
    #print r
    alpha = atang(xy_p)
    #print alpha
    beta  = alpha - ang
    beta  *=pi/180.00
    #print beta
    tras_cos = cos(beta)
    #print tras_cos
    tras_sin = sin(beta)
    n = len(xy_p[0,:])
    x = np.zeros(n)
    x += r
    x *= tras_cos
    y = np.zeros(n)
    y += r
    y *= tras_sin 

    return array([x,y])

def ss(xy_p):
    ss_x = []
    ss_y = []
    for i in xrange(len(xy_p[0,:])):
        if xy_p[1,i] >= 0.0:
            ss_x.append(xy_p[0,i])
            ss_y.append(xy_p[1,i])

    return np.array([ss_x,ss_y])


def ps(xy_p):
    ps_x = []
    ps_y = []
    for i in xrange(len(xy_p[0,:])):
        if xy_p[1,i]<= 0.0:
            ps_x.append(xy_p[0,i])
            ps_y.append(xy_p[1,i])

    return np.array([ps_x,ps_y])

def ps_ss(xy_p):
    x_le = xy_p[0,0]
    y_le = xy_p[1,0]
    ps_x = []
    ps_y = []
    ss_x = []
    ss_y = []
    ss_x.append(x_le)
    ss_y.append(y_le)
    ps_x.append(x_le)
    ps_y.append(y_le)
    
    i  = 1
    while (len(ps_x) < 2 or len(ss_x) < 2):
        
        if xy_p[1,i] >= y_le:
            ss_x.append(xy_p[0,i])
            ss_y.append(xy_p[1,i])
        else:
            ps_x.append(xy_p[0,i])
            ps_y.append(xy_p[1,i])

        i += 1
    while i < len(xy_p[0,:]) - 1:
#        y1_ss = (ss_y[len(ss_y)-1] - ss_y[len(ss_y) - 2])/(ss_x[len(ss_x)-1] - ss_x[len(ss_x) - 2])
#        y_ss_exp = ss_y[len(ss_y)-1] + y1_ss*(xy_p[0,i]- ss_x[len(ss_x) -1])
#        y1_ps = (ps_y[len(ps_y)-1] - ps_y[len(ps_y) - 2])/(ps_x[len(ps_x)-1] - ps_x[len(ps_x) - 2])
#        y_ps_exp = ps_y[len(ps_y)-1] + y1_ps*(xy_p[0,i]- ps_x[len(ps_x) -1])
 
        
        if abs(ss_y[len(ss_y) -1] - xy_p[1,i]) < abs(ps_y[len(ps_y) -1] - xy_p[1,i]):
            ss_x.append(xy_p[0,i])
            ss_y.append(xy_p[1,i])
        else:
            ps_x.append(xy_p[0,i])
            ps_y.append(xy_p[1,i])
        i += 1

    ss_x.append(xy_p[0,len(xy_p[0,:]) -1])
    ss_y.append(xy_p[1,len(xy_p[0,:]) -1])
    ps_x.append(xy_p[0,len(xy_p[0,:]) -1])
    ps_y.append(xy_p[1,len(xy_p[0,:]) -1])

        
        
    print ss_y, ps_y 
    return np.array([ss_x, ss_y]), np.array([ps_x, ps_y])

def normalize(xy_p, chord):
    
    xy_p /= chord
    return xy_p
    
def sort_ar_x (a):
    
    a = zip(*a)
    a.sort(key = lambda x : x[0])
    return np.array(zip(*a))




#if __name__=='__main__':
#    script, filename = argv
def import_prof(filename):
    xy_p = read_xy_p(filename)
    xy_te = xy_p[:,0]
    xy_le = xy_p[:,find_le(xy_p, xy_te)]
    xy_p = shift_o(xy_p, xy_le)
    ang = atang(xy_te - xy_le)
    xy_p = rotate(xy_p, ang)
    #ut.plot(xy_p[0,:], xy_p[1,:])
    chord = distance_2p(xy_te,xy_le)
    xy_p  = normalize(sort_ar_x(xy_p), chord)
    print xy_p
#    xy_ps = ps(xy_p)
#    print xy_ps
#    xy_ss = ss(xy_p)
#    print xy_ss
    xy_ss, xy_ps  = ps_ss(xy_p)
#    ut.plot(xy_ps[0,:], xy_ps[1,:])
#    ut.plot(xy_ss[0,:], xy_ss[1,:])
#    plt.show()
#    raw_input('\>')
    return xy_ss, xy_ps

#plot(x,y)
#show()


