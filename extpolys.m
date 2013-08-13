%EXTPOLYS   MCSC external (unbounded) polygonal domain.
%
%   P = extpolys(...
%     [ w_(1,1), w(2,1), ... ], ...
%     [ w_(2,1), w(2,2), ... ], ...
%     [ w_{3,1), w(3,2), ... ], ...
%   ) creates an unbounded polygonal domain where w_(k,j) represents the
%   kth vertex of the jth polygon.
%
%   P = extpolys(...
%     polygon([ w_(1,1), w(2,1), ... ]), ...
%     polygon([ w_(2,1), w(2,2), ... ]), ...
%     polygon([ w_{3,1), w(3,2), ... ]), ...
%   ) is the same as above.
%
%   See also polygons.
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

% Copyright Everett Kropf, 2011

classdef extpolys < polygons

  properties
    beta;
  end
  
  methods
    function P = extpolys(varargin)
      P = P@polygons(varargin{:});
      
      P.beta = P.vl;
      for j = 1:P.m
        P.beta(1:P.vc(j),j) = 1 - P.p{j}.alpha;
      end
    end
    
    function plot(P)
      toHold = ishold;
      hold on
      
      [m vc vl] = polydat(P);
      pbox = [inf -inf inf -inf];
      for j = 1:m
        w = [ vl(1:vc(j),j); vl(1,j) ];
        plot(real(w),imag(w),'k')
        plot(real(vl(1,j)),imag(vl(1,j)),'ro')
        pbox([1 3]) = min(pbox([1 3]),[min(real(w)) min(imag(w))]);
        pbox([2 4]) = max(pbox([2 4]),[max(real(w)) max(imag(w))]);
      end
      
      set(gca,'dataaspectratio',[1 1 1]);
      scale = [-0.7,0.7];
      pbox(1:2) = mean(pbox(1:2)) + diff(pbox(1:2))*scale;
      pbox(3:4) = mean(pbox(3:4)) + diff(pbox(3:4))*scale;
      axis(pbox)
      
      if ~toHold
        hold off
      end
    end % plot
    
  end %methods
  
end
