%FPEXTREFL   MCSC reflection method integrand for the external map
%   (unbounded) case.
%
%   I = fpextrefl(P, C, N) creates a FPEXTREFL object where P is a EXTPOLYS
%   object, C is a CIRCDOMAIN object, and N is the reflection truncation
%   parameter. Registers a listener with the circle domain object
%   (reflections are only recalculated on a change in the circle domain).
%
%   For details of the reflection algorithm, see [cite paper here].
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

% Copyright Everett Kropf 2013,
% except as may be noted below.
% Parts of this file are due to Tom DeLillo, Wichita State, 2003.

classdef fpextrefl < fprefl

  properties
    P; C; N;
    cnu; rnu; znu;
    snu; jlr; rsum;
    sing;
    
    mthdstr = 'reflection';
  end
  
  methods
    function I = fpextrefl(P,C,N)
      if nargin > 0
        I.P = P;
        I.C = C;
        I.N = N;
        
        build_reflections(I);
        
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
            b = a.eval_fkj(s(1).subs{:});
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
      build_reflections(I);
    end
    
    function fpz = eval_fprime(I,z)
      [m vc vl beta] = polydat(I.P);
      z_nu = I.znu; c_nu = I.cnu; s_nu = I.snu;
      
      n = size(s_nu,1); % numel(snu(:,1));

      zprod1 = ones(size(z));
      logzprod2 = 0;
      for j = 1:m
        for nu = 1:n
          zs = (z - c_nu(nu,j));
          for k = 1:vc(j)
            logzprod2 = logzprod2 + ...
              beta(k,j)*log(1-(z_nu(k,nu,j)-c_nu(nu,j))./zs);
          end
          zprod1=zprod1.*((zs./(z-s_nu(nu,j))).^2);
        end
      end

      fpz = zprod1.*exp(logzprod2);
      
    end % eval_fprime
    
    function fkj = eval_fkj(I,z,k,j)
      z_nu = I.znu; s_nu = I.snu; c_nu = I.cnu;
      
      n = size(s_nu,1);
      logprod = zeros(size(z));
      for nu = 1:n
        logprod = logprod ...
          + log( 1 - (z_nu(k,nu,j) - c_nu(nu,j))./(z - c_nu(nu,j)) ) ...
          + log( 1 - (c_nu(nu,j) - s_nu(nu,j))./(z - s_nu(nu,j)) );
      end
      
      fkj = exp(logprod);
    end % eval_fkj
    
    function build_reflections(I)
      m = I.P.m; vc = I.P.vc;
      
      zn = prevertices(I.C);
      
      % just mark prevertices as possible singularities for now...
      I.sing = zeros(sum(vc),1);
      b = vc(1);
      I.sing(1:b) = zn(1:b,1);
      for j = 2:m
        I.sing(b+1:b+vc(j)) = zn(1:vc(j),j);
        b = b + vc(j);
      end
      
      zn = reshape(zn,max(vc),1,m);
      cn = I.C.c; rn = I.C.r;
      cn = cn(:).'; rn = rn(:)';
            
      [zn cn rn sn jlr rsum] = I.reflectzsmi(zn,cn,rn,cn,vc,I.N);
      
      I.znu = zn; I.cnu = cn; I.rnu = rn;
      I.snu = sn; I.jlr = jlr; I.rsum = rsum;
    end
    
  end % methods
    
end
