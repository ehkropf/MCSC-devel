% interior (bounded) radial slitmap
% see [reference goes here]
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

classdef intradslitmap < intslitmap

  methods
    function S = intradslitmap(a0,c,r,N)
      if nargin == 0
        sargs = {};
      else
        rcs = 0*c; %zeros(1,numel(c));
        if nargin < 4
          N = [];
        end
        sargs = { a0, c, r, rcs, N };
      end
      S = S@intslitmap(sargs{:});
      
    end % constructor
    
  end % methods
  
end
