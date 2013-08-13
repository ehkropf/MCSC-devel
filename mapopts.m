%MAPOPTS   MCSC map option class.
%   This value class manages properties that control the behavior of
%   general MCSC map creation (things common to both bounded and unbounded
%   methods.
%     [List properties here]
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

classdef mapopts

  properties
    N;
    fpclassh;
    method;
    
    monitor = 0;
    fignum = [];
    
    tolerance = 1e-12;
    ngj = 12;
    halfrule = true;
    
    nosolver = false;
    nolinktol = false;
    
    numcontin = struct('hmin',1e-12);
    
    % KLUDGE?
    objfun = [];
  end
  
  properties (Abstract)
    % Cell array describing methods for implementing the integrand.
    % Each row is a method. The columns are (1) method string,
    % (2) function handle to integrand sublass constructor, and
    % (3) default value for the "truncation parameter" N.
    fpclasstable;
  end
  
  properties (Access=private)
    isInternal = false;
  end
    
  methods
    function o = mapopts
      % default to first fprime method in table
      o.method = o.fpclasstable{1,1};
    end % constructor
    
    function o = set.method(o,value)
      if o.isInternal, o.method = value; return, end
      n = find(strcmp(value,{o.fpclasstable{:,1}}));
      if ~any(n)
        error('mcsc:badmethod','Method %s not recognized.',value);
      end
      o.isInternal = true;
      [ o.method o.fpclassh o.N ] = deal(o.fpclasstable{n,:});
      o.isInternal = false;
    end
    
    function o = set.fpclassh(o,value)
      if o.isInternal, o.fpclassh = value; return, end
      
      vcls = class(value);
      if ~strcmp(vcls,'function_handle')
        error('mcsc:notfcnh','%s is not a valid function handle.',vcls);
      end
      
      n = [];
      vstr = func2str(value);
      for j = 1:size(o.fpclasstable,1)
        if strcmp(vstr,func2str(o.fpclasstable{j,2}))
          n = j;
        end
      end
      
      if isempty(n)
        notok = false;
        try
          tmpo = value();
        catch %#ok<CTCH>
          notok = true;
        end
        if ~notok && isa(tmpo,'integrand')
          o.fpclassh = value;
          o.isInternal = true;
          o.method = 'unknown';
          o.isInternal = false;
        else
          error('mcsc:badfprime','%s is not a valid fprime.',vstr);
        end
      else
        o.isInternal = true;
        [ o.method o.fpclassh o.N ] = deal(o.fpclasstable{n,:});
        o.isInternal = false;
      end
    end
    
    function o = set.halfrule(o,value)
      if ~islogical(value) || ~(numel(value) == 1)
        error('mcsc:badvalue','Halfrule must be true or false.')
      end
      o.halfrule = value;
    end
    
    function o = set.tolerance(o,value)
      o.tolerance = value;
      
      if ~o.isInternal && ~o.nolinktol
        o.isInternal = true;
        o.ngj = max(ceil(-log10(value)),4);
        o.isInternal = false;
      end
    end
    
    function o = set.ngj(o,value)
      o.ngj = value;
      
      if ~o.isInternal && ~o.nolinktol
        o.isInternal = true;
        o.tolerance = 10^(-value);
        o.isInternal = false;
      end
    end
    
  end
  
end
