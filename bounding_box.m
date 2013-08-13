function box = bounding_box(w,s)
% BOUNDING_BOX - get bounding box for figure
%
%   box = bounding_box(w) given a set of points  in the plane w, gives a
%      vector box = [xmin xmax ymin ymax]. Given a scale factor s = 1.2;
%      a bounding square is given by
%        s*max{ max(real(w))-min(real(w)), max(imag(w))-min(imag(w)) }
%
%   box = bounding_box(w,s) specifies s.
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

% Copyright Everett Kropf, 2013

if nargin < 2
  s = 1.2;
end

box([1 3]) = min([real(w(:)) imag(w(:))]);
box([2 4]) = max([real(w(:)) imag(w(:))]);
db = max(diff(box(1:2)),diff(box(3:4)));
box(1:2) = mean(box(1:2)) + s/2*db*[-1 1];
box(3:4) = mean(box(3:4)) + s/2*db*[-1 1];
