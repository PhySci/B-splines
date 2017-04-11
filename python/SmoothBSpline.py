"""
Class description should be here

F. Mushenok, https://github.com/PhySci 
"""

import numpy as np
from scipy import sparse
from scipy.sparse import linalg


class SmoothBSpline():

    # interpolation coefficients
    _a = np.array([0])
    _b = np.array([0])
    _c = np.array([0])
    _d = np.array([0])

    # position and weights of knots
    _x = np.array([0])
    _y = np.array([0])
    _w = np.array([0])

    # smooth parameter
    _p = 0.1

    @classmethod
    def __init__(self, x, y, w, p):
        """
        Init object and calculate spline coefficients
        :param x: array of independent variable
        :param y: array of dependent variable
        :param w: array of weights 
        :param p: smoothness of spline in range [0,1]
        """
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
        # [xi,yi,sizeval,w,origint,tol,tolred] = chckxywp(x,y,max(2,m),w,tol,'adjtol');

        # @TODO kostil'
        if True:
            self._x = x
            self._y = y
            self._w = w

        h = np.diff(self._x)
        n = h.shape[0]

        # create D (weight matrix)
        D = sparse.spdiags(self._w, 0, n+1, n+1)

        # create Q matrix
        ih = np.divide(1, h)
        print 'ih'
        print ih
        Q = sparse.spdiags([ih[:-1], -ih[1:] - ih[:-1], ih[1:]], [-1, 0, 1], n+1, n-1)

        # create  T matrix
        T = sparse.spdiags([h[1:]/3.0, 2*(h[:-1]+h[1:])/3.0, h[1:]/3.0], [-1, 0, 1], n-1, n-1)

        self._c = linalg.inv(Q.transpose() * D * D * Q + p*T) *p*Q.transpose() * self._y

        self._a = self._y - D*D*Q*self._c/self._p

        self._c = np.insert(self._c, 0, 0)
        self._c = np.append(self._c,0)
        self._d = (self._c[1:]-self._c[:-1])/(3*h)

        self._b = (self._a[1:]-self._a[:-1])/h -self._c[:-1]*h -self._d*h**2

        #print 'h=', h
        print 'b=', self._b
        print 'a=', self._a
        print 'c=', self._c
        print 'd=', self._d


    @classmethod
    def calcCurve(self, xArr):
        """
        Calculate spline curve
        :param xArr: x position of calculation points
        :return: array of y value
        """

         # sort x values
        xArr = np.sort(xArr)

        # init output array
        yArr = np.zeros(xArr.shape)

        # find bins number
        bins = np.digitize(xArr, self._x)-1
        #print 'Bins:',bins

        it = np.nditer(xArr, flags=['f_index'])

        while not it.finished:
            knotId = bins[it.index]
            diff = it[0]-self._x[knotId]
            #print xArr[it.index], self._x[knotId], diff
            yArr[it.index] = self._a[knotId] +self._b[knotId]*diff +self._c[knotId]*diff**2+self._d[knotId]*diff**3
            it.iternext()

        return yArr


    @classmethod
    def getCoeffs(self):
        """
        Return coefficients of splines
        :return: array of array
        """
        return [self._a, self._b, self._c, self._d]