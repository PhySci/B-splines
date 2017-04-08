import numpy as np
import math

# x - input x array
# y - input y array
# w - input weight array
#

def bspl(x,y,w=None,p=0.1):
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
    #[sp,values,rho] = spaps1(x,y,tol,w,m);

    if w is None:
        w = np.ones(x.shape[0],1)  # @TODO ugly piece


    # @TODO check imput arrays here
    # - size of x and y are equal
    # - remove NaNs and Inf (simulteneously in x,y and w)
    # x should be unique. If x is not unique, then process it (calculate mean x taking to account weigth arrays)
    # - anything else?
    # - p should be in [0,1]


    #
    # spaps1

    #[xi,yi,sizeval,w,origint,tol,tolred] = chckxywp(x,y,max(2,m),w,tol,'adjtol');

    # @TODO костыль на время отсуствия проверки
    #  надо потом понять как всё это парситься и приводится в порядок
    if True:
        xi = x
        yi = x
        sizeval = xi.shape[0]
        tol = 1
        tolread =1


    dx = np.diff(xi)

    """
    dx = diff(xi); n = size(xi,1); xi = reshape(xi,1,n);
    yd = size(yi,2); # yd =1
 Set up the linear system for solving for the B-spline coefficients of
  % the m-th derivative of the smoothing spline (as outlined in the note
  % [C. de Boor, Calculation of the smoothing spline with weighted roughness
  % measure] obtainable as smooth.ps at ftp.cs.wisc.edu/Approx), making use
  % of the sparsity of the system.
  %
  %  A  is the Gramian of B_{j,x,m}, j=1,...,n-m,
  % and  Ct  is the matrix whose j-th row contains the weights for
  % for the `normalized' m-th divdif (m-1)! (x_{j+m}-x_j)[x_j,...,x_{j+m}]

   % Deal with the possibility that a weighted roughness measure is to be used.
   % This is quite easy since it amounts to multiplying the integrand on the
   % interval (x(i-1) .. x(i)) by 1/tol(i), i=2:n.
   if length(tol)==1, dxol = dx;
   % elseif length(tol)==n: Actually, this is already checked in CHCKXYWP
   else
      lam = reshape(tol(2:end),n-1,1); tol = tol(1); dxol = dx./lam;
   end

# if m = 2
A = spdiags([dxol(2:n - 1), 2 * (dxol(2:n - 1)+dxol(1:n - 2)), dxol(1:n - 2)], ...
                                                                               - 1:1, n - m, n - m) / 6;
odx = 1. / dx;
Ct = spdiags([odx(1:n - 2), -(odx(2:n - 1)+odx(1:n - 2)), odx(2:n - 1)], ...
0:m, n - m, n);

% Now determine  f  as the smoothing spline, i.e., the minimizer of
  %
  %  rho*E(f) + F(D^m f)
  %
  % with the smoothing parameter  RHO  chosen so that  E(f) <= TOL .
  % Start with  RHO=0 , in which case  f  is polynomial of order  M  that
  % minimizes  E(f) .
  % If the resulting  E(f)  is too large, follow C. Reinsch
  % and determine the proper  rho  as the unique zero of the function
  %
  %   g(rho):= 1/sqrt(E(rho)) - 1/sqrt(TOL)
  %
  % (since  g  is monotone increasing and is close to linear for larger RHO)
  % using Newton's method  at  RHO = 0 but deviating from Reinsch's advice by
  % using the Secant method after that.
  % This requires  g'(rho) = -(1/2)E(rho)^{-3/2} DE(rho) , with  DE(rho)
  % derived from the determining equations for  f . These are
  %
  %        Ct y = (Ct W^{-1} C + rho A) u,  u := c/rho,
  %
  % with  c  the B-coefficients of  D^m f , in terms of which
  %
  %       y - f = W^{-1} C u,  E(rho) =  (C u)' W^{-1} C u ,
  % hence
  %         DE(rho) =  2 (C u)' W^{-1} C Du ,
  % with
  %           - A u = (Ct W^{-1} C + rho A) Du
  %
  % In particular, DE(0) = -2 u' A u , with  u = (Ct W^{-1} C) \ (Ct y) ,
  % hence  g'(0) = E(0)^{-3/2} u' A u.

   cty = Ct*yi; wic = spdiags(1./w(:),0,n,n)*(Ct.'); ctwic = Ct*wic;

must_integrate = 1;

% we are to work with a specified rho
rho = -tol;
u = (ctwic + rho * A)\cty; ymf = wic * u; values = (yi - ymf).';

if must_integrate
      sp = spmak(xi,(rho*u).');
      if exist('lam','var') % we must divide D^m s by lambda before integration.
                      % This may require additional knots at every jump
                      % in lambda.
         lam = reshape(lam,1,n-1);
         if m==1
            sp = spmak(xi,((rho*u).')./repmat(lam,yd,1));
         else
            jumps = 1+find(diff(lam)~=0); mults = [];
            if ~isempty(jumps) % must insert knots at these points
              sp = sprfn(sp, ...
                     reshape(xi(ones(m-1,1),jumps),1,(m-1)*length(jumps)));
               mults = [jumps(1)-1,diff(jumps)+m-1];
            end

          % Then must group the coefs accordingly, and divide by their LAM,
          % in anticipation of which the vector MULTS has already been
          % generated. Here are the details:
          % There are  length(jumps)+1  such groups. The first one goes
          % from 1 to jumps(1)-1, to be divided by lam(1),
          % i.e., jumps(1)-1 terms are involved. Next goes
          % from jumps(1) to jumps(2)-1+m-1, divided by lam(jumps(1)+1)
          % i.e., jumps(2)-jumps(1)+(m-1)  terms are involved. Next goes
          % from jumps(2)+m-1 to jumps(3)-1+2(m-1), by lam(jumps(2)+1)
          % i.e., jumps(3)-jumps(2)+m-1  terms are involved. Next goes
          %  ...
          % from jumps(end)+1-m+(m-1)*(end-1) to length(coefs,2), by lam(end)
          % We'll use BRK2KNT to produce the needed multiplicities.

            [knots, coefs] = spbrk(sp);
            lams = brk2knt(lam([1 jumps]), [mults,size(coefs,2)-sum(mults)]);
            sp = spmak(knots, coefs./repmat(lams,yd,1));

         end % of case distinction based on value of m
      end  % of division of D^m s by lam

      for j=1:m-1, sp = fnint(sp); end
      sp = fnint(sp,values(:,1));
   end

    % At this point, SP differs from the answer by a polynomial of order M , and
    % this polynomial is computable from its values   VALUES-FNVAL(SP,XI)
    if m>1
    [knots, coefs, ignored, k] = spbrk(sp);
    knotstar = aveknt(knots,k); knotstar([1 end]) = knots([1 end]);
    % (special treatment of endpoints to avoid singularity of collocation
    %  matrix due to noise in calculating knot averages)

    if yd==1
       vals = polyval(polyfit(xi-xi(1),values-fnval(sp,xi),m-1),knotstar-xi(1));  <-- only this part is used
    else
      %  Unfortunately, MATLAB's POLYVAL and POLYFIT only work for scalar-valued
      %  functions, hence we must make our own homegrown fit here.
      vals = fnval(spap2(1,m,xi,values-fnval(sp,xi)),knotstar);
    end

    if m==2
      sp = spmak(knots, coefs+vals); % vals give the value at the Greville
                                     % points of a straight line, and these
                                     

                                     % we know therefore to be the B-coeffs
                                     % of that straight line wrto knots.
    elseif m==3
      sp = spmak(knots, coefs+spbrk(spapi(knots,knotstar,vals),'coefs'));
    end
    end
    
    return 0
    """
    return 0