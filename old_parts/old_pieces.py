# old piece of cst.py

def class_function(x,N):
    return x**N[0]*(1-x)**N[1]

def bin_coef(i,n):
    return factorial(n)/(factorial(i)*factorial(n-i))

def cst(x,A):
    
    N = A[0:2]
    A = A[2:len(A)]
    n = len(A)
    S = 0
    for i in xrange(n):
        S += A[i]*bin_coef(i,n-1)*(1-x)**(n-1-i)*x**i

    return class_function(x,N)*S





def ps_err(ycst, yp, c = 1.0):
    r = np.where(yp != 0.0, abs(yp-ycst)/abs(yp)**c, abs(yp-ycst))
    return r

def ss_err(ycst, yp, c = 1.0):
    r = np.where(yp != 0.0, abs(yp-ycst)/abs(yp)**c, abs(yp-ycst))
    return r

def err_sum(err):
    s = 0
    for row in err:
        s +=row
    return s

def of_ss(A):

    global xy_ss
    global c_ss
    cst_ss = cst(xy_ss[0,:], A)
    err_ss = ss_err(cst_ss, xy_ss[1,:], c = c_ss)
    return  np.linalg.norm(err_ss, ord = 2.0) 
#    return err_sum(err_ss)

def of_ps(A):

    global xy_ps
    global c_ps
    cst_ps = cst(xy_ps[0,:], A)
    err_ps = ps_err(cst_ps, xy_ps[1,:], c = c_ps)
    return  np.linalg.norm(err_ps, ord = 2.0) 
#    return err_sum(err_ps)
