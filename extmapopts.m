%EXTMAPOPTS   MCSC external (unbounded) map options.
%   Subclasses MAPOPTS, specifically defines the fpclasstable property to
%   control how f' is calculated (integrand selection). See
%   [cite papers here] for details on the different methods.
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

% Copyright Everett Kropf, 2011, 2012

classdef extmapopts < mapopts

  properties
    % See comment in MAPOPTS for adding to this table
    fpclasstable = {
      'refl' @fpextrefl 6
      'fpls' @fpextfpls 25
    };
  end

  methods
    function o = extmapopts
      o = o@mapopts;
    end % constructor
    
  end
  
end
