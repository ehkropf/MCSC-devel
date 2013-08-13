%OBJECTIVEFUN   MCSC objective function base class.
%  Base class for creating objective functions (bounded or unbounded).
%  Encapsulates some things common to all MCSC objective functions such
%  as calculating the separation parameter, and common monitor and plot
%  funtionality. OBJECTIVEFUN is a handle subclass.
%
%  O = objectivfun(P, monitor, fignum) creates an OBJECTIVEFUN object where
%  P is a polygons container object, monitor is (treated as) a boolean
%  value indicating solution monitoring will occur, and fignum is the
%  figure number the objective function plots to.
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

classdef objectivefun < handle

  properties
    separam = 0;
    seplim = 0;
    
    fh;
    ah;
    ntimes = 0;
    monitor = 0;
  end
  
  properties (Abstract)
    P; % polygons object
    C; % circdomain object
  end

  methods
    function O = objectivefun(P, monitor, fignum)
      if nargin > 0 && ~isempty(P)
        O.seplim = 1/nthroot(P.m-1,4);
        
        if nargin > 1 && ~isempty(monitor)
          O.monitor = monitor;
          if monitor
            if nargin > 2 && ~isempty(fignum)
              O.fh = figure(fignum);
            else
              O.fh = figure;
            end
            
            clf(O.fh);
            O.ah = gca;
          end
        end
        
      end
    end % constructor
    
    function calc_separam(O,c,r)
      [Ri Rj] = meshgrid(r);
      [Ci Cj] = meshgrid(c);
      
      mu = (Ri+Rj)./abs(Ci-Cj);
      mu(mu==inf) = 0;
      O.separam = max(mu(:));
    end
    
    function do_monitor(O)
      O.ntimes = O.ntimes + 1;
      plot(O)
    end
    
    function plot(O)
      if O.monitor
        set(0,'currentfigure',O.fh);
        set(O.fh,'currentaxes',O.ah);
        cla
      end
    end
    
  end
end
