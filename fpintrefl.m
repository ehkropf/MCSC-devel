%FPINTREFL   MCSC reflection method integrand for the internal map
%   (bounded) case.
%
%   I = fpintrefl(P, C, N) creates a FPINTREFL object where P is a INTPOLYS
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

% Copyright Everett Kropf, 2013,
% except as may be noted below.
% Parts of this file are due to Tom DeLillo, Wichita State, 2003.

classdef fpintrefl < fprefl

  properties
    P; C; N;
    cnu; rnu; znu;
    snu; jlr; rsum;
    sing;
    
    mthdstr = 'reflection';
  end
  
  methods
    function I = fpintrefl(P,C,N)
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
      m = I.P.m; vc = I.P.vc; beta = I.P.beta;
      zn = I.znu; cn = I.cnu;
      N = I.N;
      
      zlogsum = zeros(size(z));
      
      for k = 1:vc(1)
        zlogsum = zlogsum + beta(k,1)*log(1 - z/zn(k,1,1));
      end
      nub = 0;
      nu1 = 1;
      for level=1:N
        nua = nub+1;
        if m ~= 2
          nub = ((m-1)^level - 1)/(m-2);
        elseif m == 2
          nub = nua;
        end
        for j=2:m
          for nu = nua:nub
            nu1 = nu1 + 1;
            for k=1:vc(j)
              zlogsum = zlogsum + beta(k,j)*log(1 - ...
                (zn(k,nu,j) - cn(nu,j))./(z - cn(nu,j)));
            end
            
            for k=1:vc(1)
              zlogsum = zlogsum + beta(k,1)*log(1 - ...
                (zn(k,nu1,1) - cn(nu,j))./(z - cn(nu,j)));
            end
          end
        end
      end
      
      fpz = exp(zlogsum);
    end % eval_fprime
    
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
      
      % reflect interior circles to exterior
      for j = 2:m
        [cn(j) rn(j)] = I.centrad(cn(1),rn(1), cn(j),rn(j));
        zn(1:vc(j),1,j) = I.reflect(cn(1),rn(1), zn(1:vc(j),1,j));
      end
      
      [zn cn rn sn jlr rsum] = I.reflectzsmi(zn,cn,rn,cn,vc,I.N);
      
      I.znu = zn; I.cnu = cn; I.rnu = rn;
      I.snu = sn; I.jlr = jlr; I.rsum = rsum;
    end
    
  end % methods
  
end
