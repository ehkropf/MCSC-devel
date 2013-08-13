%INTMAPOPTS   MCSC internal (bounded) map options.
%   Subclasses MAPOPTS, specifically defines the fpclasstable property to
%   control how f' is calculated (integrand selection). The bounded case
%   currently only supports the reflection method [cite paper here].
%
%   See also maptops.
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

classdef intmapopts < mapopts

  properties
    % See comment in MAPOPTS for adding to this table
    fpclasstable = {
      'refl' @fpintrefl 6
    };
  end

  methods
    function o = intmapopts
      o = o@mapopts;
    end % constructor
    
  end
  
end
