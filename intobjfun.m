%INTOBJFUN   MCSC internal (bounded) objective function.
%   OBJECTIVEFUN (and handle) subclass implementing the objective function
%   for the bounded map parameter problem described in
%   [cite paper by author here].
%
%   O = extobjfun(P, C) creates an EXTOBJFUN object where P is a polygons
%   container object and C is a circdomain object.
%
%   O = extobjfun(P, C, opts) reads a mapopts object to get monitor and
%   fignum information to pass to the parent constructor.
%
%   See also objectivefun, mapopts.
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

classdef intobjfun < objectivefun

  properties
    P;
    C;
    N;
    z0;
    w0;
    quad;
  end
  
  methods
    function O = intobjfun(P,fixp,C,opts)
      if nargin > 0 && ~isempty(P)
        args = { P };
        
        if nargin > 3 && ~isempty(opts)
          args = { args{:} opts.monitor, opts.fignum };
        end
      else
        args = {};
      end
      O = O@objectivefun(args{:});
      
      if nargin > 0
        O.P = P;
        O.C = C;
        O.N = opts.N; %N;
        O.z0 = fixp(1);
        O.w0 = fixp(2);
%         O.seplim = 1/nthroot(O.P.m-1,4);
        
        % fprime determined by class handle in options object
        fpo = opts.fpclassh(P,C,opts.N);
        O.quad = gjquad(P,fpo,opts.ngj,opts.halfrule);
      end
    end % constructor
    
    function b = subsref(a,s)
      % impliments following behavior:
      %   F = a(X) -- as function
      %   a.property
      %   a.property.subproperty
      %   a.property.subproperty(index)
      
      errcond = 0;
      
      switch s(1).type
        case '()'
          % zp = a(z) -- treat as function
          b = a.fun_eval(s(1).subs{:});
        case '.'
          if length(s) > 1
            if length(s) > 2
              switch s(3).type
                case '()'
                  % a.prop.subprop(stuff)
                  b = a.(s(1).subs).(s(2).subs)(s(3).subs{:});
                otherwise
                  errcond = 1;
              end
            else
              % a.prop.subprop
              b = a.(s(1).subs).(s(2).subs);
            end
          else
            % a.prop
            b = a.(s(1).subs);
          end
        otherwise
          errcond = 1;
      end
      
      if errcond
        error('MCSC:invaliduse',...
          'Invalid use of object of class %s.',class(a));
      end
    end % subsref
    
    function F = fun_eval(O,Xu)
      % circle domain does constr/unconstr x-form internally
      O.C.Xu = Xu;
      
      if O.monitor
        do_monitor(O);
      end

      qobj = O.quad;
      [m,vc,vl] = polydat(O.P);
      c = O.C.c; r = O.C.r;
      tkj = O.C.t;
      zkj = prevertices(O.C);
      
      F = zeros(sum(vc) + 3*m-4, 1);
            
      % scaling/rotation constant
      Q12 = arcq(qobj, tkj(1,1), 1, tkj(2,1), 2, 0, 1);
      A = (vl(2,1) - vl(1,1))/Q12;
      
      % side lengths
      SL = zeros(max(vc),m);

      for j = 1:m
        off = sum(vc(1:j-1));
        if j < 3
          % don't need the last sidelength on C1 or C2
          lpvn = off + (1:vc(j)-1)'; rpvn = lpvn + 1;
          left = tkj(1:vc(j)-1,j); right = [left(2:end); tkj(vc(j),j)];
        else
          lpvn = off + (1:vc(j))'; rpvn = [lpvn(2:end); off+1];
          left = tkj(1:vc(j),j); right = [left(2:end); tkj(1,j)];
        end        
        right(left>right) = right(left>right) + 2*pi;
        cv = repmat(c(j),size(lpvn,1),1);
        rv = repmat(r(j),size(lpvn,1),1);
        
        SL(1:size(lpvn,1),j) = ...
          A*arcq(qobj, left,lpvn, right,rpvn, cv,rv);
      end
      
      
      % paths from z0 to circle prevertices z_(1,j)
      left = zkj(1,1:m).'; right = O.z0*ones(m,1);
      pvn = 1 + cumsum([0 vc(1:j-1)])';
      T = -A*lineq(qobj, left,pvn, right,0*pvn);
      
      
      % orientation + sidelength conditions of first sides for circles
      % 2,...,m
      for j = 2:m
        F(2*j-3) = real( SL(1,j) - (vl(2,j) - vl(1,j)) );
        F(2*j-2) = imag( SL(1,j) - (vl(2,j) - vl(1,j)) );
      end
      
      % position conditions for circles 1,...,m wrt w0
      for j = 1:m
        F(2*(m-1)+2*j-1) = real( T(j) - (vl(1,j) - O.w0) );
        F(2*(m-1)+2*j) = imag( T(j) - (vl(1,j) - O.w0) );
      end
      
      % sidelength conditions
      jk = 4*m-2;
      
      % leave off last side lengths of circles 1 and 2
      for j = 1:2
        for k = 2:vc(j)-1
          F(jk+k-1) = abs(SL(k,j)) - abs(vl(k+1,j) - vl(k,j));
        end
        jk = jk + vc(j)-2;
      end
      
      % all sidelengths of remaining
      for j = 3:m
        for k = 2:vc(j)-1
          F(jk+k-1) = abs(SL(k,j)) - abs(vl(k+1,j) - vl(k,j));
        end
        F(jk+vc(j)-1) = abs(SL(vc(j),j)) - abs(vl(1,j) - vl(vc(j),j));
        jk = jk + vc(j)-1;
      end
      
    end % fun_eval
    
    function calc_separam(O)
      c = O.C.c; r = O.C.r;
      for j = 2:O.P.m
        [c(j) r(j)] = fprefl.centrad(c(1),r(1),c(j),r(j));
      end
      
      calc_separam@objectivefun(O,c,r);
    end
    
    function plot(O)
      toHold = ishold;
      
      plot@objectivefun(O)
      
      plot(O.C)
      hold on
      
      m = O.P.m;
      zkj = prevertices(O.C);
      plot(real(O.z0),imag(O.z0),'rx')
      for j = 1:m
        lin = [O.z0; zkj(1,j)];
        plot(real(lin),imag(lin),'k-')
      end
      
      calc_separam(O);
      
      title(sprintf(['objective function call #%d\n'...
                     '(m-1)^{-1/4} = %5.4f     \\Delta = %5.4f'],...
                     O.ntimes,O.seplim,O.separam))
      
      drawnow
      
      if ~toHold
        hold off
      end
    end
    
  end % methods
  
end
