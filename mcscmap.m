%MCSCMAP   MCSC map superclass.
%  This class is incomplete and must be subclassed. It encapsulates some of
%  the things common to both bounded and unbounded maps:
%    - creating a quadrature object
%    - high-level map evaluation operations including
%      - finding nearest circle to a point
%      - calculating the integration multiplicative constant
%  Each subclass needs in its constructor to create an objective function
%  and a solver to be used here. In addition subclasses must define an
%  integration path method to be used in map evaluation. MCSCMAP is a value
%  class.
%
%  M = mcscmap(P, N, C, O, S) is the superclass constructor. P is the
%  polygon container object, N is the truncation parameter (various
%  interpretations, see below), C is the circdomain container object, O is
%  an objectivefun object, and S is a solver object. This constructer
%  gives some informative output, creates a quadrature object, and calls
%  the S.run_solver method.
%
%  The meaning of the truncation parameter will change based on the
%  numerical representation of the map. See papers by author.
%
%  mcscmap methods:
%    integpath   (abstract method) integration path used in map_eval;
%                required to be defined by subclasses.
%    map_eval    evaluates the map at a given point
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

classdef mcscmap

  properties
    P; % polygonal domain
    C; % circle domain
    N; % truncation parameter
    A; % multiplicative constant
    fp; % MCSC integrand (fprime)
    Q; % MCSC quadrature object
    S; % solver
  end
  
  methods
    function M = mcscmap(P,N,C,O,S,varargin)
      if nargin > 0
        % !this should check the credentials of its arguments
        % !do we need varargin?
        
        M.Q = O.quad;
        M.fp = M.Q.fp;

        cstr = upper(class(M));
        fprintf(['\n========================================'...
          '================================================\n']);
        
        if isempty(S)
          fprintf('%s: creating map using the %s method with N=%d.\n\n',...
            cstr,M.fp.mthdstr,N);
        else
          fprintf(['%s: solving the parameter problem using '...
            'the %s method with N=%d.\n'],cstr,M.fp.mthdstr,N)
          fprintf('%s: tolerance = %.0e, ngj = %d\n\n',...
            cstr,S.tolerance,M.Q.ngj);
        
          run_solver(S);
        end
                        
        M.P = P;
        M.C = C;
        M.N = N;
        M.S = S;
      end
    end % constructor
    
    function b = subsref(a,s)
      % implements the following behavior:
      %   w = a(z) -- as a function (map_eval)
      %   a.property
      
      errcond = 0;
      
      switch s(1).type
        case '()'
          % treat as function
          b = a.map_eval(s(1).subs{:});
        case '.'
%             b = a.(s(1).subs);
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
          errcond = 1;
      end
      
      if errcond
        error('MCSC:invaliduse',...
          'Invalid use of object of class %s.',class(a));
      end
    end % subsref
    
    function w = map_eval(M,z)
      A = M.A; qobj = M.Q; vl = M.P.vl;
      c = M.C.c; r = M.C.r;
      
      w = nan(size(z));
      for v = 1:numel(z)
        if isnan(z(v)), continue, end
        
        [arc,line,k,j] = integpath(M,z(v));
        
        w(v) = vl(k,j);
        
        if ~isempty(arc)
          w(v) = w(v) + A*arcq(qobj, arc.left,arc.lpvn, ...
            arc.right,0, c(j),r(j));
        end
        if ~isempty(line)
          w(v) = w(v) + A*lineq(qobj, line.left,line.lpvn, z(v),0);
        end
      end
    end % map_eval

    function [j,t] = nearestcircle(M,z)
      % find the nearest circle to a given point z
      % basically from what Driscoll had in 2005
      
      c = M.C.c; r = M.C.r;
      
      d = abs(z - c);
      [dmin,j] = min(abs(d - r));
      t = mod(angle(z - c(j)),2*pi);
    end
        
    function A = calc_constant(M)
      qobj = M.Q;
      vl = M.P.vl;
      tkj = M.C.t;
      
      A = (vl(2,1) - vl(1,1)) / arcq(qobj, tkj(1,1),1, tkj(2,1),2, 0,1);
    end % calc_constant
        
  end % methods
  
  methods (Abstract)
    [arc,line,k,j] = integpath(M,z);
  end % abstract methods
  
end
