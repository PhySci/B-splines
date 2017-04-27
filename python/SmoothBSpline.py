"""
Class description should be here

F. Mushenok, https://github.com/PhySci 
"""

import numpy as np
import matplotlib.pyplot as plt

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

    _knots = np.array([]) # array of knots
    _coeffs = np.array([]) # array of coefficients
    _k = 0 # order of spline

    @classmethod
    def __init__(self, x=0, y=0, w=0, p=0):
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
    def bspl(self, x, y, w, p):
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
        #x = np.array(x)
        #y = np.array(y)

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
        tol = p # np.array([1])
        tolread = 1

        dx = np.diff(xi)
        n = xi.shape[0]
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

        A = sparse.spdiags(np.array([dxol[1:n-2]/6, 2*(dxol[1:n-2]+dxol[:n-3])/6, dxol[:n-3]/6]), np.array([-1, 0, 1]), n - m, n - m)

        #plt.plot(dx)
        #plt.show()
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
        self._knots = knots
        self._coeffs = coeffs
        self._k = 4

        knotstar = self.aveknt(self._knots,self._k)
        knotstar[0] = knots[0]
        knotstar[-1] = knots[-1]

        pf = np.polyfit(xi - xi[0], values - self.eval(xi), 1)
        vals = np.polyval(pf, knotstar-xi[0])

        self._coeffs = self._coeffs+vals

        return 0

    @classmethod
    def fnint(self, t, a, n, k, val = 0.0):
        """
        Integration of BB spline
        :param t: array of knots
        :param a: array of coefficients
        :param n: number of coefficients
        :param k: order of spline
        :param val: value of the integral on the left end
        :return: list of new knots and coefficients  
        """

        # increase multiplicity of last knot to k
        index =  np.flatnonzero(np.diff(t))
        needed = index[-1]+1 - n
        if (needed > 0):
            t = np.hstack([t, np.tile(t[n+k-1], needed)])
            a =np.hstack([a, np.zeros(1, needed)])
            n = n+needed

        if (val!=0.0):
            needed = k - index[0]-1
            knots = np.append(np.hstack((np.tile(t[0],[1, needed+1]).ravel(), t)), t[n+k-1])
            slice = np.hstack((val,np.zeros(needed), np.multiply(a,np.tile(np.divide(t[k:k+n]-t[:n],k), 1))))
            coeffs = np.cumsum(slice)
        else:
            knots= np.hstack((t, t[n+k-1]))
            coeffs =  np.cumsum(np.multiply(a, np.tile( np.divide(t[k:k+n] - t[:n], k), 1)))

        return knots, coeffs

    @classmethod
    def eval(self,x,left = 0):
        """
        :param x: array of absciss where function will be calculated 
        :param left: uknown parameter @TODO sort out what it is
        :return: array of function values for each x value
        """
        nx = 1 #everywhere because 1D array of x values is accepted
        mx = x.shape[0]
        lx = mx * nx;
        xs = x # reshape(x, 1, lx);

        # Take apart spline:
        t = self._knots
        a = self._coeffs
        k = self._k
        d = 1
        n = a.shape[0]

        # If there are no points to evaluate at return empty matrix of appropriate size:
        if (lx == 0):
            v = np.zeros(0)
            return v
        #Otherwise, augment the knot sequence so that first and last knot each has multiplicity >= K.
        #(AUGKNT would not be suitable for this since any change in T must be accompanied by a corresponding change in A.)

        # I think this peace of code does not work. Keep it it comments
        #index = np.flatnonzero(np.diff(t))
        #addl = k - index[0]
        #addr = index[-1] - n
        #if (addl > 0 || addr > 0):
        #    npk = n + k
        #    t1 = t[np.ones(addl)]
        #    t2 = t[:npk-1]
        #    t3 = t[npk]
        #    t4 =
        #    t = t[, 0:npk-1, npk[ones(1, addr)]]
        #    a = [zeros(d, addl) a zeros(d, addr)]
        #    n = n + addl + addr

        index = np.digitize(x, t[k-1:n])  # [~, index] = histc(xs, [-inf, t(k + 1:n), inf]);
        NaNx = np.flatnonzero(index == 0)
        index = np.minimum(index+k-1,n) #quite important line, debug it carefully

        if NaNx.shape[0]>0:
            index[NaNx] = k

        # Now, all is ready for the evaluation.
        # Carry out in lockstep the first spline evaluation algorithm
        # (this requires the following initialization):
        if (k > 1):
            dindex = index
            tx = t[np.tile(np.arange(2-k,k), [lx,1]) + np.tile(dindex,[2*(k-1),1]).transpose() -1]
            tx = tx - np.tile(xs, [2*(k - 1),1]).transpose()
            b = np.tile(np.arange(d*(1-k),1,d), [lx, 1]) + np.tile(index, [k, 1]).transpose()
            b = a[b-1];

        # (the following loop is taken from SPRPP)
        for r in range(1,k):
            for i in range(1,k-r):
                a5 = np.multiply(tx[:, i+k-2],b[:, i-1]) - np.multiply(tx[:, i+r-2],b[:, i])
                a6 = tx[:,i+k-2] - tx[:,i+r-2]
                b[:, i-1] = np.divide(a5,a6)

        v = b[:,0]
        #Finally, zero out all values for points outside the basic interval:
        v[x<t[0]] = 0
        v[x>t[-1]] = 0
        return v


    @classmethod
    def aveknt(self,t, k):
        """
        AVEKNT Knot averages.
        Returns the averages of successive K-1 knots, i.e. the points
    
              TSTAR(i) = (T_{i + 1} + ... + T_{i + K - 1} ) / (K - 1)
    
        Recommended as good interpolation point choices when interpolating from  S_{K, T}.
        For example, with k and the increasing sequence  breaks  given, the statements
        
        t = augknt(breaks, k);
        x = aveknt(t);
        sp = spapi(t, x, sin(x));
    
        provide spline interpolant to the sine function on the interval [breaks(1)..breaks(end)].
    
        :param t: 
        :param k: 
        :return: 
        """
        t = self._knots;
        n = t.shape[0] - k;
        if (k < 2):
            print 'SPLINES:AVEKNT:wrongk'
            return 0
        elif (n < 0):
            print 'SPLINES:AVEKNT:toofewknots'
            return 0
        elif (k == 2):
            tstar = np.reshape(t[np.arange(1,n+1)], [1, n])
        else:
            temp = np.tile(t, [k-1, 1]).transpose().flatten(1)

            c = np.hstack([temp,np.zeros(k-1)])
            c1 = np.reshape(c, [k-1, n+k+1])
            c2 = c1.sum(0)
            temp = np.divide(c2,k-1)
            tstar = temp[1:n+1]
        return tstar

    def checkData(self,x,y,w,p = 0.1):
        """
        Check and adjust input
        :param x: 
        :param y: 
        :param nmin: 
        :param w: 
        :param p: 
        :param adjtol: 
        :return: 


        """


        # make sure X is real
        #if ~all(isreal(x))
        #    x = real(x);
        #    warning(message('SPLINES:CHCKXYWP:Xnotreal'))
        #end

        # deal with NaN's and Inf's among the sites:
        #    nanx = find(~isfinite(x));
        #if ~isempty(nanx)
        #    x(nanx) = [];
        #    warning(message('SPLINES:CHCKXYWP:NaNs'))
        #end

        #n = length(x);

        # re - sort, if needed, to ensure nondecreasing site sequence:
        tosort = False;
        if np.any(np.diff(x) < 0):
            tosort = True
            ind = np.argsort(x)
            x = x[ind]
            y = y[ind]
            w = w[ind]

        #nstart = n + length(nanx);


        # make sure that sites, values and weights match in number:


        if (x.shape!=y.shape):
            raise ValueError("X don't match Y")

        if (x.shape!=w.shape):
            raise ValueError("W don't match Y")


        #Remove values and error weights corresponding to nonfinite sites:
        nanf = np.logical_or(np.logical_or(~np.isfinite(x),~np.isfinite(y)), ~np.isfinite(w))
        if np.any(nanf):
            x = x[~nanf]
            y = y[~nanf]
            w = w[~nanf]



        """ 
        if ~isempty(nanx), y(nanx,:) = []; if nonemptyw, w(nanx) = [];
        end
        if roughnessw % as a first approximation, simply ignore the
        % specified weight
        to the left of any ignored point.
        p(max(nanx, 2)) = [];
        end
        end
        if tosort, y = y(ind,:); if
        nonemptyw, w = w(ind);
        end, end

        % deal
        with nonfinites among the values:
            nany = find(sum(~isfinite(y), 2));
        if ~isempty(nany)
            y(nany,:) = [];
            x(nany) = [];
            if nonemptyw, w(nany) =[]; end
            warning(message('SPLINES:CHCKXYWP:NaNs'))
            n = length(x);
            if n < minn
                error(message('SPLINES:CHCKXYWP:toofewX', sprintf('%g', minn))), end
            if roughnessw % as a first approximation, simply ignore the
            % specified
            weight
            to
            the
            left
            of
            any
            ignored
            point.
        p(max(nany, 2)) = [];
        end
        end

        if nargin == 3 & & nmin, return, end % for SPAPI, skip the averaging

        if nargin > 3 & & isempty(w) % use the trapezoidal rule weights:
            dx = diff(x);
            if any(dx), w = ([dx;0]+[0;dx]).'/2;
            else, w = ones(1, n);
            end
            nonemptyw = ~nonemptyw;
        end

        tolred = 0;
        if ~all(diff(x)) % conflate repeat sites, averaging the corresponding values
        % and summing
        the
        corresponding
        weights
        mults = knt2mlt(x);
        for j=find(diff([mults;0]) < 0).'
        if nonemptyw
            temp = sum(w(j - mults(j):j));
            if nargin > 5
                tolred = tolred + w(j - mults(j):j)*sum(y(j - mults(j):j,:).^ 2, 2);
                end
                y(j - mults(j),:) = (w(j - mults(j):j)*y(j - mults(j):j,:)) / temp;
                w(j - mults(j)) = temp;
                if nargin > 5
                    tolred = tolred - temp * sum(y(j - mults(j),:).^ 2);
                    end
                else
                    y(j - mults(j),:) = mean(y(j - mults(j):j,:), 1);
                    end
                end

                repeats = find(mults);
                x(repeats) = [];
                y(repeats,:) = [];
                if nonemptyw, w(repeats) =[]; end
                if roughnessw % as a first approximation, simply ignore the
                % specified
                weight
                to
                the
                left
                of
                any
                ignored
                point.
            p(max(repeats, 2)) = [];
            end
            n = length(x);
            if n < minn, error(message('SPLINES:CHCKXYWP:toofewX', sprintf( '%g', minn ))), end
        end

        if nargin < 4, return, end

        % remove
        all
        points
        corresponding
        to
        relatively
        small
        weights(since
        a
        % (near -)
        zero
        weight in effect
        asks
        for the corresponding datum to be dis-
        % regarded
        while , at the same time, leading to bad condition and even
        % division
        by
        zero).
        origint = []; % this
        will
        be
        set
        to
        x([1 end]).
        ' in case the weight for an end
        % data
        point is near
        zero, hence
        the
        approximation is computed
                         % without
        that
        endpoint.
        if nonemptyw
        ignorep = find(w <= (1e-13) * max(abs(w)));
        if ~isempty(ignorep)
           if ignorep(1) == 1 | | ignorep(end) == n, origint = x([1 end]).
        '; end
        x(ignorep) = [];
        y(ignorep,:) = [];
        w(ignorep) = [];
        if roughnessw
           % as a
        first
        approximation, simply
        ignore
        the
        % specified
        weight
        to
        the
        left
        of
        any
        ignored
        point. \
            p(max(ignorep, 2)) = [];
        end
        n = length(x);
        if n < minn
        error(message('SPLINES:CHCKXYWP:toofewposW', sprintf('%g', minn)))
        end
        end
        end
"""
        return x, y, w