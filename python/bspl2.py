import numpy as np
from scipy import sparse
from scipy.sparse import linalg

import matplotlib.pylab as plt
import math


# x - input x array
# y - input y array
# w - input weight array
#

def bspl2(x, y, w=None, p=0.1, verbose = False):
    """
    Return smooth B spline
    :param x: array of independent variable
    :param y: array of dependent variable
    :param w: array of weights 
    :param p: smoothness of spline in range [0,1]
    :param verbose: detailed output of function
    :return: still don't know @TODO
    """
    m = 2

    # matlab code
    # order of curve (I suppose it means cubic)
    # [sp,values,rho] = spaps1(x,y,tol,w,m);

    # Convert to numpy array
    x = np.array(x)
    y = np.array(y)

    if w is None:
        # @TODO ugly piece
        w = np.ones([x.shape[0], 1])

    # @TODO check imput arrays here
    # - size of x and y are equal
    # - remove NaNs and Inf (simulteneously in x,y and w)
    # x should be unique. If x is not unique, then process it (calculate mean x taking to account weigth arrays)
    # - anything else?
    # - p should be in [0,1]


    #
    # spaps1

    # [xi,yi,sizeval,w,origint,tol,tolred] = chckxywp(x,y,max(2,m),w,tol,'adjtol');

    # @TODO kostil'
    if True:
        xi = x
        yi = y


    h = np.diff(xi)
    n = h.shape[0]

    # create D (weight matrix)
    D = sparse.spdiags(w,0,n+1,n+1)

    # create Q matrix
    ih = np.divide(1,h)
    Q = sparse.spdiags([ih[:-1],-ih[1:]-ih[:-1],ih[:-1]],[-1,0,1],n+1,n-1)

    # create  T matrix
    T = sparse.spdiags([h[:-1]/3.0,2.0*(h[:-1]+h[1:])/3.0,h[:-1]/3.0],[-1,0,1],n-1,n-1)

    c = linalg.inv(Q.transpose()*D*D*Q + p*T)*p*Q.transpose()*y
    a = y - (1.0/p)*D*D*Q*c

    c = np.insert(c,0,0)
    c[-1] = 0
    d = (c[1:]-c[:-1])/h[:-1]
    b = (a[1:-1]-a[:-2])/h[:-1]  - c[-1]*h[-1] - d*h[:-1]*h[:-1]

    print 'h=',h
    print 'b=', b
    print 'a=', a
    print 'c=', c
    print 'd=', d

    it = np.nditer(x,flags=['f_index'])

    y2 = np.zeros(n+1)
    for ind in it:
        y2[ind] = a[ind]*1.0 #+b[ind]+c[ind]+d[ind]

    return y2;

if __name__ == '__main__':

    x = [0, 1, 2, 3, 4, 5]
    y =  [0, 1, 8, 27, 64, 125]
    print y
    w = [1, 1/3, 1/5, 1/15, 10, 1]

    k = bspl2(x,y,w,0.1)
    print k

    plt.plot(x, y,'rx')
    plt.plot(x, k, '-go')

    plt.show()
