%FPEXTFPLS   MCSC finite product/least squares method integrand for the
%   external (unbounded map.
%
%   I = fpextfpls(P, C, N) creates a FPEXTFPLS object. See FPEXTREFL for
%   the arguments.
%
%   For details of the algorithm, see [cite paper here].
%
%   See also fpextrefl.
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

% Copyright Everett Kropf, 2011

classdef fpextfpls < integrand

  properties
    P; C; N;
    M; x; sing;
    
    side_cond;
    
    mthdstr = 'least squares';
  end
  
  properties (Access=private)
    c;
    r;
    zkj;
  end
  
  methods
    function I = fpextfpls(P,C,N)
      if nargin > 0
        I.P = P;
        I.C = C;
        I.N = N;
        
        % makes the system to find coefficients just overdetermined
        I.M = 2*N+1;
        
        I.c = C.c; I.r = C.r; I.zkj = prevertices(C);
        calc_coefficients(I);
        
        % add listener for circdomain update via Xu (done from the
        % objective function, other updates of the circle domain not
        % supported -- call build_reflections manually)
        addlistener(C,'externXuUpdate',@I.handleCircDomChg);
      end
    end % constructor
    
    function b = subsref(a,s)
      % impliments following behavior:
      %   zp = a(z) -- as function
      %   a.property
      %   a.property.subproperty
      %   a.property.subproperty(index)
      
      errcond = false;
      
      switch s(1).type
        case '()'
          % zp = a(z) -- treat as function
          if numel(s(1).subs) > 1
            b = a.eval_xkj(s(1).subs{:});
          else
            b = a.eval_fprime(s(1).subs{:});
          end
        case '.'
          if length(s) > 1
            if length(s) > 2
              switch s(3).type
                case '()'
                  % a.prop.subprop(stuff)
                  b = a.(s(1).subs).(s(2).subs)(s(3).subs{:});
                otherwise
                  errcond = true;
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
         errcond = true;
      end
      
      if errcond
        error('MCSC:invaliduse',...
          'Invalid use of object of class %s.',class(a));
      end
    end % subsref
    
    function handleCircDomChg(I,src,evnt) %#ok<INUSD>
      I.c = I.C.c; I.r = I.C.r; I.zkj = prevertices(I.C);
      calc_coefficients(I);
    end
    
    function fpz = eval_fprime(I,z)
      [m vc vl beta] = polydat(I.P);
      
      fpz = ones(size(z));
%       logfpz = zeros(size(z));
      for j = 1:m
        for k = 1:vc(j)
          fpz = fpz .* ( eval_xkj(I,z,k,j) ).^beta(k,j);
%           logfpz = logfpz + beta(k,j)*eval_xkj(I,z,k,j);
        end
      end
%       fpz = exp(logfpz);
    end % eval_fprime
    
    function calc_coefficients(I)
      [m vc] = polydat(I.P);
      
      x = zeros(m*I.N,sum(vc));
      for j = 1:I.P.m
        off = sum(vc(1:j-1));
        for k = 1:vc(j)
          x(:,k+off) = make_xkj(I,k,j);
        end
      end
      I.x = x;
      
      % mark prevertices as possible singularities...
      I.sing = zeros(sum(vc),1);
      b = vc(1);
      I.sing(1:b) = I.zkj(1:b,1);
      for j = 2:m
        I.sing(b+1:b+vc(j)) = I.zkj(1:vc(j),j);
        b = b + vc(j);
      end
      
      % update side condition
      calc_side_cond(I);
    end % calc coefficients
    
    function calc_side_cond(I)
      [m vc vl beta] = polydat(I.P);
      c = I.c(:); r = I.r(:); zkj = I.zkj;
      x = I.x; N = I.N;
      
      sc = 0;
      for j = 1:m
        dP = x( N*(0:m-1)+1, (1:vc(j))+sum(vc(1:j-1)) ); %#ok<PROP>
        dP = dP.';
        sc = sc + beta(1:vc(j),j)'*( zkj(1:vc(j),j) - c(j) - dP*r ); %#ok<PROP>
      end
      
      I.side_cond = abs(sc);
    end % side condition
    
    function x = make_xkj(I,k,j)
      m = I.P.m;
      c = I.c; r = I.r;
      N = I.N; M = I.M;
      
      aj = I.zkj(k,j);
      sj = c(j);
      
      %% need M equally spaced points around each circle
      % adjust points to work around logrithmic branch cut
%       z = repmat(exp(1i*(2*pi*(1:M)'/M-pi)),1,m);
      z = repmat(exp(1i*(2*pi*(1:M)'/M)),1,m);
      z = z*diag(r);
      z = repmat(c(:).',M,1) + z;
      z = z(:);
      
      
      %% setup B1 and B2
      A = Amatrix(I,z);
      B1 = imag(A); B2 = real(A);
      
      jrows = (j-1)*M+1:j*M;
      zG = z(jrows);
      G = zeros(M,m*N);
      for p = 1:m
        G(:,(p-1)*N+1:p*N) = ...
          A(jrows,(p-1)*N+1:p*N) .* (-(zG-sj)./(zG-c(p))*(1:N));
      end
      GR = real(G); GI = imag(G);
      
      B1(jrows,:) = GR;
      B2(jrows,:) = -GI;
      
      
      %% the E matrix
      D = toeplitz([-1 zeros(1,M-2)],[-1 1 zeros(1,M-2)]); %%#ok<NASGU>
%       arg = repmat(' ',1,4+2*m);
      arg = cell(1,m);
%       sp = 1;
%       for p = 1:m
%         if p == j
%           arg(sp:sp+5) = 'eye(M)'; sp = sp+6;
%         else
%           arg(sp) = 'D'; sp = sp+1;
%         end
%         if p < m, arg(sp) = ','; sp = sp+1; end
%       end
%       E = eval(['blkdiag(' arg ')']);
      for p = 1:m
        if p == j
          arg{p} = eye(M);
        else
          arg{p} = D;
        end
      end
      E = sparse(blkdiag(arg{:}));
      
      %% right-hand side
      argas = (z-aj)./(z-sj);
      argas(jrows) = 0;
      argas = angle(argas);
      rhs = -E*argas;
            
      %% solve linear system
      B = E*[B1 B2];
      x = B\rhs;
      
%       fprintf('\nResidual (2) for k=%d, j=%d: %0.15e\n',k,j,norm(B*x-rhs));
%       fprintf('Residual (inf) for k=%d, j=%d: %0.15e\n',k,j,norm(B*x-rhs,inf));

      x = complex(x(1:end/2),x(end/2+1:end));

    end % make xkj
    
    function w = eval_xkj(I,z,k,j)
      A = Amatrix(I,z(:));
      w = zeros(size(z));
      x = I.x(:,k+sum(I.P.vc(1:j-1)));
%       w(:) = log(1 - (I.zkj(k,j)-I.c(j))./(z(:)-I.c(j))) + A*x;
%       w(:) = (1 - (I.zkj(k,j)-I.c(j))./(z(:)-I.c(j))) .* exp(A*x);
      w(:) = (z(:)-I.zkj(k,j))./(z(:)-I.c(j)) .* exp(A*x);
    end % eval xkj
    
    function A = Amatrix(I,z)
      % return matrix of basis functions for g
      
      m = I.P.m; N = I.N;
      c = I.c; r = I.r;
      
      A = zeros(length(z),m*N);
      for n = 1:m
        zcn = r(n)./(z-c(n)); %#ok<PROP>
        A(:,(n-1)*N+1) = zcn;
        for ell = 2:N
          A(:,(n-1)*N+ell) = A(:,(n-1)*N+ell-1) .* zcn;
        end
      end
    end % A matrix
    
  end % methods
    
end
