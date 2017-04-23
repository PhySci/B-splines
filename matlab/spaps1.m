function [sp,values,rho] = spaps1(x,y,tol,w,m)
%SPAPS1 Univariate smoothing spline.

[xi,yi,sizeval,w,origint,tol,tolred] = chckxywp(x,y,max(2,m),w,tol,'adjtol');

if tol(1)>=0 % we may have to reduce the tolerance, because of data points
             % with the same site but different values ...
   tol(1) = tol(1)-tolred;
   if tol(1)<0
      warning(message('SPLINES:SPAPS:toltoolow'));
      tol(1) = 0;
   end
end   

dx = diff(xi); n = size(xi,1); xi = reshape(xi,1,n);
yd = size(yi,2);

if n==m % the smoothing spline is the interpolating polynomial
   rho = 0;
   values = yi.'; sp = spmak(augknt(xi([1 n]),m),yi(ones(m,1),:).');
elseif tol(1)==0||tol(1)==(-inf)  % return the interpolating variational spline;
   rho = inf; values = yi.';
   switch m
   case 1
      sp = spmak(augknt(xi,2),yi.');
   case 2
      sp = fn2fm(csape(xi,yi.','variational'),'B-');
   case 3
      % here's a quick fix on getting the quintic variational interpolating
      % spline by getting the complete quintic spline interpolant to vector-
      % valued properly augmented data, ...
      y = yi.'; yo = zeros(4,n+4); yo(:,[2 3 end-1 end]) = eye(4);
      k = 6;  s = spapi(augknt(xi,k),xi([1 1 1:n n n]),...
         [y(:,1) zeros(yd,2) y(:,2:n) zeros(yd,2); yo]);
      % then want the spline  sp = [eye(yd) A]*s  with  A  chosen so that
      % D^j sp (x(i)) = 0  for  j=3:4, i = [1 end] , i.e.,
      %   0 = [eye(yd) A]*[s4 s3] =: [eye(yd) A]*[B;C], hence  A = -B/C :
      s3 = fnder(s,3);
      BC = [fnval(fnder(s3),xi([1 end])), fnval(s3,xi([1 end]))];
      sp = fncmb(s,[eye(yd) -BC(1:yd,:)/BC(yd+1:end,:)]);
   otherwise
      error(message('SPLINES:SPAPS:dontdosmoothingorder', num2str( m )))
   end
else

  % Set up the linear system for solving for the B-spline coefficients of
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

   if m==1
      A = spdiags(dxol,0,n-1,n-1);
      Ct = spdiags([-ones(n-1,1) ones(n-1,1)], 0:1, n-1,n);
   elseif m==2
      A = spdiags([dxol(2:n-1), 2*(dxol(2:n-1)+dxol(1:n-2)), dxol(1:n-2)],...
                                         -1:1, n-m,n-m)/6;
      odx = 1./dx;
      Ct = spdiags([odx(1:n-2), -(odx(2:n-1)+odx(1:n-2)), odx(2:n-1)], ...
                                                0:m, n-m,n);
   elseif m==3
      % have some fun putting together the Gramian by Gauss quadrature
      % in preparation for the case of general  m .
      pdx = reshape((.774596669241483/2)*dx,1,n-1);
      % get interval midpoints ...
      temp = reshape((xi(1:n-1)+xi(2:n))/2,1,n-1);
      % ... and the m Gauss points in each interval ...
      temp = [temp-pdx;temp;temp+pdx];
      % ... and, from this, the values of all the B-splines at all these points
      [ig,ig,ig,ig,blocks] = bkbrk(spcol(augknt(xi,m),m,temp(:).',1));
      % BLOCKS=[A_1; ...; A_{n-1}] consists of  n-1  m-by-m  matrices  A_i
      % with  A_i(p,q) the value at the p-th point in interval (xi(i)..xi(i+1))
      % of the q-th B-spline relevant there. Hence, ...
      blocks = reshape(blocks,m,m*(n-1));
      % ... reorganizes this into [B_1, ..., B_m], with  B_q(p,i) the value
      % at the p-th point in interval  i  of the q-th B-spline relevant there.
      % For each interval  i , we need to compute the properly weighted scalar
      % product of all B-splines relevant there, the weights being the
      % quadrature weights, for the interval of length  dx , which are
      % (5/18, 4/9, 5/18)*dx
      % Hence, multiply the p-th row of BLOCKS by the p-th weight:
      wblocks = [(5/18)*blocks(1,:);(4/9)*blocks(2,:);(5/18)*blocks(3,:)];
      % then get the 1-row matrix [C_11, C_12, C_13, C_22, C_23, C_33] with
      %  C_qr(i) the integral over the interval  i  of the (pointwise) product
      % of the relevant B-splines  q  and  r  there:
      nm1 = 1:n-1; nm3 = 1:(n-3);
      temp = reshape(...
         sum(blocks(:,[nm1 nm1 nm1 n-1+nm1 n-1+nm1 2*(n-1)+nm1]).* ...
         wblocks(:,[ 1:3*(n-1)    n:3*(n-1)     2*(n-1)+nm1])).* ...
         reshape(repmat(dxol,1,6),1,6*(n-1)),6*(n-1),1);
      % Note that SPDIAGS aligns columns (rather than rows) in this case
      % and that, e.g., the i-th main-diagonal entry is the sum of three terms,
      % namely C_33(i) + C_22(i+1) + C_11(i+2); etc
        A = spdiags([temp(2*(n-1)+2+nm3), ...
                   temp(4*(n-1)+1+nm3)+temp(n+1+nm3), ...
                   temp(5*(n-1)+nm3)+temp(3*(n-1)+1+nm3)+temp(2+nm3), ...
                   temp(4*(n-1)+nm3)+temp(n+nm3), ...
                   temp(2*(n-1)+nm3)],-2:2,n-m,n-m);

      dxx = dx(1:(n-2))+dx(2:(n-1)); dxxx = dxx(1:(n-3))+dx(3:(n-1));
      odx = 1./dx; odxodx = odx(1:(n-2)).*odx(2:(n-1));
      odxx = 1./dxx;
      % note that SPDIAGS aligns rows (rather than columns) in this case
      Ct = 2*spdiags([-odx(1:(n-3)).*odxx(1:(n-3)), ...
                    dxxx.*odxodx(1:(n-3)).*odxx(2:(n-2)), ...
                   -dxxx.*odxx(1:(n-3)).*odxodx(2:(n-2)), ...
                    odxx(2:(n-2)).*odx(3:(n-1))], 0:m, n-m,n);
   else
      error(message('SPLINES:SPAPS:cantdomoothingorder', num2str( m )))
   end

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

   cty = Ct*yi;
   wic = spdiags(1./w(:),0,n,n)*(Ct.');
   ctwic = Ct*wic;

   must_integrate = 1;
   if tol<0      % we are to work with a specified rho
      rho = -tol;
      a1 = ctwic + rho*A
      a2 = cty
      u = a1\a2;
      ymf = wic*u;
      values = (yi - ymf).';
   else          % determine  rho  from the tolerance requirement

      u = ctwic\cty; ymf = wic*u; E = trace(u'*Ct*ymf);
      if E<tol
         values = (yi - ymf).'; rho = 0;
         sp = spmak(augknt(xi([1 n]),m), values(:,ones(1,m)));
         must_integrate = 0;
      else
         oost = 1/sqrt(tol); g0 = 1/sqrt(E) - oost;
         rho = -g0*E*sqrt(E)/trace(u'*A*u); delrho = rho;
         while ~isnan(rho)&&delrho>0
            u = (ctwic + rho*A)\cty; ymf = wic*u; E = trace(u'*Ct*ymf);
            if 100*abs(E-tol)<tol, break, end
            grho = 1/sqrt(E) - oost;
            delrho = delrho/(g0/grho-1);
            g0 = grho; rho = rho+delrho;
         end
         values = (yi - ymf).';
      end
   end
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

      for j=1:m-1
          sp = fnint(sp);
      end
      sp = fnint(sp,values(:,1));
   end
end

 % At this point, SP differs from the answer by a polynomial of order M , and
 % this polynomial is computable from its values   VALUES-FNVAL(SP,XI)
if m>1
   [knots, coefs, ignored, k] = spbrk(sp);
   knotstar = aveknt(knots,k); knotstar([1 end]) = knots([1 end]);
   % (special treatment of endpoints to avoid singularity of collocation
   %  matrix due to noise in calculating knot averages)
   if yd==1
      vals = polyval(polyfit(xi-xi(1),values-fnval(sp,xi),m-1),knotstar-xi(1));
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

if ~isempty(origint), sp = fn2fm(fnbrk(fn2fm(sp,'pp'),origint),'B-'); end
if length(sizeval)>1, sp = fnchg(sp,'dz',sizeval); end