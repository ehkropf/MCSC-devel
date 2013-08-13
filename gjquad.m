% MCSC G-J quadrature
%
%
% This file is part of MCSC.
% 
% MCSC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MCSC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MCSC.  If not, see <http://www.gnu.org/licenses/>.

% Copyright Everett Kropf, 2013,
% except as may be noted below.

classdef gjquad < quadrature

  properties
    fp;
    P;
    ngj;
    nodes;
    wgts;
    betav;
    
    halfrule;
  end
  
  methods
    function Q = gjquad(P,integ,ngj,halfrule)
      if nargin > 0
        if nargin < 4 || halfrule
          Q.halfrule = true;
        else
          Q.halfrule = false;
        end
        
        if nargin < 3
          Q.ngj = 12;
        else
          Q.ngj = max(ngj,4);
        end
        
        Q.fp = integ;
        Q.P = P;
        calc_qdata(Q);
      end
    end
    
    function calc_qdata(Q)
      np = Q.ngj;
      m = Q.P.m; vc = Q.P.vc; beta = Q.P.beta;
      
      nodes = zeros(np,sum(vc)+1);
      wgts = nodes;
      betav = zeros(sum(vc),1);
      for j = 1:m
        v = sum(vc(1:j-1));
        for k = 1:vc(j)
          betav(v+k) = beta(k,j);
          if beta(k,j) <= -1, continue, end
          [nodes(:,v+k) wgts(:,v+k)] = gaussj(Q,0,beta(k,j));
        end
      end
      [nodes(:,end) wgts(:,end)] = gaussj(Q,0,0);
      Q.nodes = nodes;
      Q.wgts = wgts;
      Q.betav = betav;
    end
    
    function [x,w] = gaussj(Q,alf,bet)
      %GAUSSJ Nodes and weights for Gauss-Jacobi integration.
      %   [X,W] = GAUSSJ(N,ALF,BET) returns nodes and weights for Gauss-Jacobi
      %   integration. Z and W are N-vectors such that
      %
      %          _ +1
      %         /
      %         |               ALF     BET
      %         |     f(x) (1-x)   (1+x)      dx
      %         |
      %        _/
      %           -1
      %
      %   is approximated by sum(f(Z) .* W).
      %
      %   Copyright 1997 by Toby Driscoll. Last updated 04/11/97.
      
      %   Uses the Lanczos iteration connection to orthogonal polynomials.
      %   Borrows heavily from GAUSSJ out of SCPACK Fortran.
      
      % Calculate coeffs a,b of Lanczos recurrence relation (closed form is
      % known).  Break out n=1 specially to avoid possible divide by zero.
      n = Q.ngj;
      apb = alf+bet;
      a(1) = (bet-alf)/(apb+2);
      b(1) = sqrt(4*(1+alf)*(1+bet) / ((apb+3)*(apb+2)^2));
      N = 2:n;
      a(N) = (apb)*(bet-alf) ./ ((apb+2*N).*(apb+2*N-2));
      N = 2:(n-1);
      b(N) = sqrt(4*N.*(N+alf).*(N+bet).*(N+apb) ./ ...
        (((apb+2*N).^2-1).*(apb+2*N).^2));
      
      % Find eigvals/eigvecs of tridiag "Ritz" matrix
      if n > 1
        [V,D] = eig(diag(a) + diag(b,1) + diag(b,-1));
      else
        V = 1;
        D = a;
      end
      
      % Compute normalization (integral of w(x))
      c = 2^(apb+1)*gamma(alf+1)*gamma(bet+1)/gamma(apb+2);
      
      % return the values
      x = diag(D);
      w = c*(V(1,:)').^2;
      [x,ind] = sort(x);
      w = w(ind);
    end
    
    function value = arcq(Q,left,lpn,right,rpn,c,r)
      % left and right can be column vectors of prevertex thetas.
      % lpn and rpn give the ordered number for the prevertex, zero if not
      % a prevertex. Prevertex order numbers can be found by calling
      % gjquad.pvnum.
      % c and r are vectors of centers and radii of the arc for each
      % left-right pair.
      
      % sing is a vector of possible singularities to use with the 1/2 rule
      sing = Q.fp.sing;
      
      nI = numel(left);
      
%       % THIS SHOULD BE DONE OUTSIDE HERE!!!!
%       while any(right(:) < left(:))
%         wrong = right < left;
%         right(wrong) = right(wrong) + 2*pi;
%       end
      
      % any intervals need splitting? we integrate from a (possible)
      % singularity on the "left" to a non-singularity on the "right".
      le = left(:); re = right(:);
      lpn = lpn(:); rpn = rpn(:);
      c = c(:); r = r(:);
      split = (rpn ~= 0);
      if any(split)
        mid = (le(split) + re(split))/2;
        le = [ le; re(split) ];
        re(split) = mid;
        re = [ re; mid ];
        c = [ c; c(split) ];
        r = [ r; r(split) ];
        lpn = [ lpn; rpn(split) ];
      end
      
      
      val = zeros(size(le));
      for v = 1:length(le)
        a = le(v);
        dir = sign(re(v)-le(v));
        sng = sing(abs(sing - c(v)-r(v)*exp(1i*a)) > 100*r(v)*eps);
        
        tlen = r(v)*abs(le(v)-re(v));
        subdiv = -1;
        while r(v)*abs(le(v)-a) + 10*eps(tlen) < tlen
          subdiv = subdiv + 1;
          
          if Q.halfrule
            % 1/2 rule: how far do we step?
            len = min([r(v)*abs(a - re(v)); ...
              2*abs(c(v) + r(v)*exp(1i*a) - sng(:))]);
            b = a + dir*len/r(v);
          else
            b = re(v);
          end
          
          if a == le(v) && lpn(v)
            % do g-j integration; need to know what nodes and wgts to use!
            x = Q.nodes(:,lpn(v));
            w = Q.wgts(:,lpn(v));
            xline = (b - a)*(x+1)/2 + a;
            zarc = c(v) + r(v)*exp(1i*xline);
            farc = Q.fp(zarc);
            val(v) = val(v) + w' * ...
              ((0.5i*(b-a) * farc .* (zarc - c(v))) ...
              ./ ((1+x).^Q.betav(lpn(v))));
          else
            % regular gaussian integration
            x = Q.nodes(:,end);
            w = Q.wgts(:,end);
            xline = (b - a)*(x+1)/2 + a;
            zarc = c(v) + r(v)*exp(1i*xline);
            farc = Q.fp(zarc);
            val(v) = val(v) + w' * ...
              (0.5i*(b-a) * farc .* (zarc - c(v)));
          end
          a = b;
          
          if subdiv > 100 % arbitrary!
            error('mcsc:manysubdiv',['Too many subdivisions detected. '...
              'Left pvtx num: %d; right pvtx num: %d'],lpn(v),rpn(v));
          end

        end
      end
      
      % put stuff back together
      if any(split)
        vsre = val(nI+1:end);
        val = val(1:nI);
        val(split) = val(split) - vsre;
      end
      
      value = left;
      value(:) = val(:);
    end % arcq
    
    function value = lineq(Q,left,lpn,right,rpn)
      % left and right can be column vectors of prevertex thetas.
      % lpn and rpn give the ordered number for the prevertex, zero if not
      % a prevertex. Prevertex order numbers can be found by calling
      % gjquad.pvnum.
      
      % sing is a vector of possible singularities to use with the 
      % 1/2 rule.
      sing = Q.fp.sing;
      
      nI = numel(left);
      
      % any intervals need splitting? we integrate from a (possible)
      % singularity on the left to a non-singularity on the right.
      le = left(:); re = right(:);
      lpn = lpn(:); rpn = rpn(:);
      split = (rpn ~= 0);
      if any(split)
        mid = (le(split) + re(split))/2;
        le = [ le; re(split) ];
        re(split) = mid;
        re = [ re; mid ];
        lpn = [ lpn(:); rpn(split) ];
      end
      
      
      val = zeros(size(le));
      for v = 1:length(le)
        a = le(v);
        dir = (re(v)-le(v))/abs(re(v)-le(v));
        sng = sing(abs(sing - a) > 100*eps);
        
        tlen = abs(le(v)-re(v));
        subdiv = -1;
        while abs(le(v)-a) + 100*eps(tlen) < tlen
          subdiv = subdiv + 1;
          
          if Q.halfrule
            % 1/2 rule: how far do we step?
            len = min([abs(a - re(v)); 2*abs(a - sng(:))]);
            b = a + dir*len;
          else
            b = re(v);
          end

          if a == le(v) && lpn(v)
            % do g-j integration; need to know what nodes and wgts to use!
            x = Q.nodes(:,lpn(v));
            w = Q.wgts(:,lpn(v));
            zline = (b - a)*(x+1)/2 + a;
            fline = Q.fp(zline);
            val(v) = val(v) + w' * ...
              (0.5*(b - a)*fline ./ ((1+x).^Q.betav(lpn(v))));
          else
            % regular gaussian integration
            x = Q.nodes(:,end);
            w = Q.wgts(:,end);
            zline = (b - a)*(x+1)/2 + a;
            fline = Q.fp(zline);
            val(v) = val(v) + w' * (0.5*(b - a)*fline);
          end
          a = b;
          
          if subdiv > 100
            error('mcsc:manysubdiv',['Too many subdivisions detected. '...
              'Left pvtx num: %d; right pvtx num: %d'],lpn(v),rpn(v));
          end
        end
      end
      
      % put stuff back together
      if any(split)
        vsre = val(nI+1:end);
        val = val(1:nI);
        val(split) = val(split) - vsre;
      end
      
      value = left;
      value(:) = val(:);
    end % lineq
    
    function num = pvnum(Q,k,j)
      vc = Q.P.vc;
      for v = 1:numel(k)
        k(v) = k(v) + sum(vc(1:j(v)-1));
      end
      num = k;
    end

  end % methods
   
end
