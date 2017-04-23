"""
Class description should be here

F. Mushenok, https://github.com/PhySci 
"""

import numpy as np

from scipy import sparse
from scipy.sparse import linalg
#from scipy import linalg

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
        ""
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
        """


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


    @classmethod
    def bspl(self, x, y, w=None, p=0.1):
        """
        Return smooth B spline
        :param x: array of independent variable
        :param y: array of dependent variable
        :param w: array of weights
        :param p: smoothness of spline in range [0,1]
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
        #if True:
        xi = x
        yi = y
        sizeval = xi.shape[0]
        tol = 1 # np.array([1])
        tolread = 1

        dx = np.diff(xi)
        n = xi.shape[0]
        xi = xi.reshape(1, n)  # I don't understand why xi is reshaped
        yd = 1

        """
        Set up the linear system for solving for the B-spline coefficients of the m-th derivative of the smoothing spline
        (as outlined in the note [C. de Boor, Calculation of the smoothing spline with weighted roughness measure]
        obtainable as smooth.ps at ftp.cs.wisc.edu/Approx), making use of the sparsity of the system.


        Deal with the possibility that a weighted roughness measure is to be used. This is quite easy since it amounts to
        multiplying the integrand on the interval (x(i-1) .. x(i)) by 1/tol(i), i=2:n.
        if length(tol)==1, dxol = dx;
        else
          lam = reshape(tol(2:end),n-1,1); tol = tol(1); dxol = dx./lam;
        end
        """
        #if tol.shape[0] == 1:
        dxol = dx
        #else: # I suppose that part of code will never work
        #    lam = tol[2:].reshape(n - 1, 1)
        #    tol = tol[0]
        #    dxol = np.divide(dx, lam)

        """"
        A  is the Gramian of B_{j,x,m}, j=1,...,n-m,. It is an almost diagonal square matrix with (n-m) side and columns
        [dxol(2:n - 1), 2 * (dxol(2:n - 1)+dxol(1:n - 2)), dxol(1:n - 2)] on diagonals (-1, 0, 1)
        """

        A = sparse.spdiags(np.array([dxol[1:n-2]/6, 2*(dxol[1:n-2]+dxol[:n-3])/6, dxol[:n-3]/6]), np.array([-1, 0, 1]), n - m, n - m) #MATLAB TESTED

        odx = np.divide(1, dx)
        """
        Ct  is the matrix whose j-th row contains the weights for for the `normalized'
        m-th divdif (m-1)! (x_{j+m}-x_j)[x_j,...,x_{j+m}]
        
        Ct = spdiags([odx(1:n-2), -(odx(2:n-1)+odx(1:n-2)), odx(2:n-1)], 0:m, n-m,n);
        
        """

        # Ct matrix is wrong
        Ct = sparse.diags([odx[:n-2], -odx[:n-2]-odx[1:n-1], odx[1:n-1]], [0, 1, 2], shape = [n-2, n])

        """
        Now determine  f  as the smoothing spline, i.e., the minimizer of
            rho*E(f) + F(D^m f)
        with the smoothing parameter  RHO  chosen so that  E(f) <= TOL. Start with  RHO=0 , in which case  f  is polynomial
        of order  M  that minimizes  E(f). If the resulting  E(f)  is too large, follow C. Reinsch and determine the proper
        rho as the unique zero of the function
            g(rho):= 1/sqrt(E(rho)) - 1/sqrt(TOL)
        (since  g  is monotone increasing and is close to linear for larger RHO) using Newton's method  at  RHO = 0
         but deviating from Reinsch's advice by using the Secant method after that. This requires
            g'(rho) = -(1/2)E(rho)^{-3/2} DE(rho)
         with  DE(rho) derived from the determining equations for  f . These are
            Ct y = (Ct W^{-1} C + rho A) u,  u := c/rho
         with  c  the B-coefficients of  D^m f , in terms of which
             y - f = W^{-1} C u,  E(rho) =  (C u)' W^{-1} C u ,
         hence DE(rho) =  2 (C u)' W^{-1} C Du, with  - A u = (Ct W^{-1} C + rho A) Du

         In particular, DE(0) = -2 u' A u , with  u = (Ct W^{-1} C) \ (Ct y), hence  g'(0) = E(0)^{-3/2} u' A u.
         """

        cty = Ct*yi

        wic = np.dot(sparse.diags(np.divide(1, w), 0, shape = [n, n]), Ct.transpose())

        ctwic = np.dot(Ct,wic)
        # ctwic looks nice

        must_integrate = 1;

        # we are to work with a specified rho
        # @TODO fix it
        rho = tol

        # values of u vector is slightly diffre from MatLab values. Perhaps, due to another solving method.
        u = linalg.spsolve(ctwic+rho*A,cty)
        ymf = wic*u
        values = (yi - ymf).transpose()
        c = rho *u
        c = np.transpose(c)

        [knots, coeffs] = self.fnint(xi, c, c.shape[0], 2)
        [knots, coeffs] = self.fnint(knots, coeffs, coeffs.shape[0], 3, values[0])
        print 'knots', knots
        print 'coeffs', coeffs


        """ 
        % At this point, SP differs from the answer by a polynomial of order M , and
        % this polynomial is computable from its values   VALUES-FNVAL(SP,XI)
        if m>1
        """
        k = 3
        # knotstar = aveknt(knots,k) @TODO
        #knotstar[1 end] = knots[1 end]

        """
        %(special treatment of endpoints to avoid singularity of collocation
        % matrix due to noise in calculating knot averages)
        
        polyfit(x,y,n) finds the coefficients of a polynomial p(x) of degree n that fits the data, p(x(i)) to y(i),
        in a least squares sense.
        """

        a1 = xi - xi(1) # смещение на ноль
        a2 = values - fnval(sp, xi) # похоже на функцию ошибки

        # нужен вот этот кусок
        # --->> fnval(sp, xi)

         a3 = knotstar-xi(1) # смещение на ноль

        vals = polyval(np.polyfit(a1,a2,1), a3)

        """
        % vals give the value at the Greville points of a straight line, and these
        % we know therefore to be the B-coeffs of that straight line wrto knots.
        sp = spmak(knots, coefs+vals); 

        return 0
        """
        return 0

    @classmethod
    def fnint(self, t, a, n, k, val = 0.0):
        # increase multiplicity of last knot to k
        index =  np.flatnonzero(np.diff(t))

        needed = index[-1]+1 - n
        if (needed > 0):
            t = np.hstack([t, np.tile(t[:,n+k-1], [1, needed])])
            a =np.hstack([a, np.zeros(1, needed)])
            n = n+needed

        if (val!=0.0):
            needed = k - index[0]-1
            knots = np.append(np.hstack((np.tile(t[0],[1, needed+1]).ravel(), t)), t[n+k-1])
            slice = np.hstack((val,np.zeros(needed), np.multiply(a,np.tile(np.divide(t[k:k+n]-t[:n],k), 1))))
            coeffs = np.cumsum(slice)
        else:
            knots= np.hstack((t[0,:], t[0,n+k-1]))
            coeffs =  np.cumsum(np.multiply(a, np.tile( np.divide(t[0,k:k+n] - t[0,:n], k), 1)))

        return knots, coeffs