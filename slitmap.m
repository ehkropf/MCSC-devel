% radial/circular slit map superclass
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

classdef slitmap

  properties
    newtol = 1e-6;
    maxiter = 20;
  end

  methods
    function b = subsref(a,s)
      % implements the following behavior:
      %   w = a(z) -- as a function (map_eval)
      %   a.property
      %   a.property(arg)
      
      errcond = 0;
      
      switch s(1).type
        case '()'
          % treat as function
          b = a.map_eval(s(1).subs{:});
        case '.'
          if numel(s) == 1
            b = a.(s(1).subs);
          else
            b = a.(s(1).subs)(s(2).subs{:});
          end
        otherwise
          errcond = 1;
      end
      
      if errcond
        error('MCSC:invaliduse',...
          'Invalid use of object of class %s.',class(a));
      end
    end % subsref
    
    function z = mapinv_odeig(S,w,z0)
      [n,m] = size(w);
      w = w(:); z0 = z0(:);
      
      scale = w - S.map_eval(z0);
      y0 = [real(z0); imag(z0)];
      [t,y] = ode45(@(t,y) invodefun(t,y), [0 0.5 1], y0);
      
      z = complex(y(end,1:numel(w)), y(end,numel(w)+1:end));
      z = reshape(z,n,m);
      
      
      function dy = invodefun(t,y) %#ok<INUSL>
        z = complex(y(1:length(y)/2), y(length(y)/2+1:end));
        dz = scale./S.prime(z);
        dy = [real(dz); imag(dz)];
      end
      
    end % mapinv odeig
    
    function z = mapinv_newt(S,w,z0)
      [n,m] = size(z0);
      zn = z0(:);
      w = w(:);
      
      done = false(size(zn));
      done(isnan(zn)) = true;
      
      k=0;
      while ~all(done) && k < S.maxiter
        F = S.map_eval(zn(~done)) - w(~done);
        zn(~done) = zn(~done) - F./S.prime(zn(~done));
        done(~done) = abs(F) < S.newtol;
        k = k+1;
      end
      
      z = reshape(zn,n,m);
    end % mapinv newt
    
  end % methods
  
end
