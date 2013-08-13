% multipy connected bounded radial/circular slit map
% see [need reference here]
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

classdef intslitmap < slitmap

  properties
    a0;
    c;
    r;
    rcs;
    N = 30;
    x;
  end
  
  methods
    function S = intslitmap(a0,c,r,rcs,N)
      if nargin > 0
        S.a0 = a0;
        S.c = c(:).';
        S.r = r(:)';
        S.rcs = rcs;
        
        if nargin > 4 && ~isempty(N) && N
          S.N = N;
        end
        
        S = calc_coefficients(S);
      end
    end % constructor
    
    function S = calc_coefficients(S)
      a0 = S.a0; c = S.c; r = S.r; N = S.N;
      rcs = S.rcs;
      m = numel(c);

      %% sample points
      M = 2*N + 1;
      z = repmat(c,M,1) + repmat(exp(2i*pi*(0:M-1)'/M),1,m)*diag(r);
      z = z(:);
      
      %% find coefficients
      F = zeros(m*M,2*m*N);
      Fp = zeros(M,m*N);
      rh = zeros(m*M,1);
      za0 = z - a0;
      
      for p = 1:m
        blk = (p-1)*M+1:p*M;
        
        for j = 1:m
          if j == 1
            zcj = z(blk);
            Fp(:,1) = zcj;
          else
            zcj = r(j)./(z(blk)-c(j)); %#ok<PROP>
            Fp(:,(j-1)*N+1) = zcj;
          end
          for ell = 2:N
            Fp(:,(j-1)*N+ell) = Fp(:,(j-1)*N+ell-1) .* zcj;
          end
        end
        
        if p == 1 || rcs(p) %#ok<PROP>
          % circular slit
          F(blk,:) = [ real(Fp) -imag(Fp) ];
          rh(blk) = -log(abs(za0(blk)));
        else
          % radial slit
          F(blk,:) = [ imag(Fp) real(Fp) ];
          rh(blk) = -angle(za0(blk));
          if any(diff(rh(blk)) > pi)
            rh(blk) = mod(rh(blk),2*pi);
          end
        end
      end
      
      E = sparse(kron(eye(m),toeplitz([-1 zeros(1,M-2)],[-1 1 zeros(1,M-2)])));
      
      x = (E*F)\(E*rh);
      
      S.x = complex(x(1:numel(x)/2),x(numel(x)/2+1:end)); %#ok<PROP>
    end % calc coefficients
    
    function A = Amatrix(S,z)
      N = S.N; c = S.c; r = S.r;
      m = numel(c);
      
      A = zeros(numel(z),m*N);
      for jj = 1:m
        if jj == 1
          zcj = z(:);
          A(:,1) = zcj;
        else
          zcj = r(jj)./(z(:)-c(jj)); %#ok<PROP>
          A(:,(jj-1)*N+1) = zcj;
        end
        for nn = 2:N
          A(:,(jj-1)*N+nn) = A(:,(jj-1)*N+nn-1) .* zcj;
        end
      end
    end % Amatrix
    
    function w = map_eval(S,z)
      A = Amatrix(S,z);
      w = zeros(size(z));
      w(:) = (z(:) - S.a0).*exp(A*S.x);
    end % map eval
    
    function wp = prime(S,z)
      N = S.N; c = S.c; r = S.r;
      m = numel(c);
      
      A = Amatrix(S,z);
      B = zeros(numel(z),m*N);
      for jj = 1:m
        if jj == 1
          zcj = z(:);
          B(:,1) = 1;
        else
          zcj = r(jj)./(z(:)-c(jj)); %#ok<PROP>
          B(:,(jj-1)*N+1) = -zcj.*(zcj/r(jj)); %#ok<PROP>
        end
        for nn = 2:N
          B(:,(jj-1)*N+nn) = B(:,(jj-1)*N+nn-1)*nn/(nn-1) .* zcj;
        end
      end
      
      wp = zeros(size(z));
      wp(:) = exp(A*S.x) .* (1 + (z(:)-S.a0) .* (B*S.x));
    end % prime
    
    function targ = pretip_angles(S)
      % calculate pre-tip angles on circles by minimizing real and imaginary
      % parts of log(f(z)) for z = c+r*exp(1i*theta) wrt theta

      % some sample points
      m = numel(S.c);
      npts = 100;
      theta = 2*pi*(0:npts-1)'/(npts-1);
      ww = repmat(S.c,npts,1) + repmat(exp(1i*theta),1,m)*diag(S.r);
      ww = map_eval(S,ww);
      
      % common part of the derivative
      innerdt = @(z,j) (z-S.c(j)).*prime(S,z)./map_eval(S,z);

      targ = zeros(2,m-1);
      for j = 2:m
        if S.rcs(j)
          % circular slits
          mmw = angle(ww(:,j));
          if any(abs(diff(mmw))>pi)
            mmw = mod(mmw,2*pi);
          end
          
          dt = @(t) real(innerdt( S.c(j)+S.r(j)*exp(1i*t), j));
        else
          % radial slits
          mmw = abs(ww(:,j));
          
          dt = @(t) -imag(innerdt( S.c(j)+S.r(j)*exp(1i*t), j));
        end
        
        [t adx] = min(mmw);
        [t bdx] = max(mmw);
        targ(:,j-1) = [ fzero(dt,theta(adx)); fzero(dt,theta(bdx)) ];
      end
      
    end % pretip_angles
    
  end % methods
  
end
