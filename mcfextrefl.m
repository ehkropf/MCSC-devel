% MCSC factor, exterior case by reflections
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
% Parts of this file are due to Tom DeLillo, WSU, 2003.

classdef mcfextrefl < handle

  properties
    C; N;
    m; vc;
  end
  
  methods
    function F = mcfextrefl(P,C,N)
      if nargin > 0
        F.m = P.m;
        F.vc = P.vc;
        F.C = circdomain(C); % copy circle domain
        F.N = N;
      end
      
      % build reflections
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
          b = a.eval_factor(s(1).subs{:});
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
    
    function fpz = eval_factor(F,z)
      [m vc vl beta] = polydat(I.P);
      znu = I.znu; cnu = I.cnu; snu = I.snu;
      
      n = numel(snu(:,1));

      zprod1 = ones(size(z));
      logzprod2 = 0;
      for j = 1:m
        for nu = 1:n
          zs = (z - cnu(nu,j));
          for k = 1:vc(j)
            logzprod2 = logzprod2 + ...
              beta(k,j)*log(1-(znu(k,nu,j)-cnu(nu,j))./zs);
          end
          zprod1=zprod1.*((zs./(z-snu(nu,j))).^2);
        end
      end

      fpz = zprod1.*exp(logzprod2);
      
    end % eval_factor
    
  end

end
