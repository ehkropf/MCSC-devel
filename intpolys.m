%INTPOLYS   MCSC interior (bounded) polygonal domain.
%
%   P = intpolys(...
%     [ w_(1,1), w(2,1), ... ], ...
%     [ w_(2,1), w(2,2), ... ], ...
%     [ w_{3,1), w(3,2), ... ], ...
%   ) creates a bounded polygonal domain where w_(k,j) represents the
%   kth vertex of the jth polygon. The first polygon must enclose the
%   others.
%
%   P = intpolys(...
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

% Copyright Everett Kropf, 2013,
% except as may be noted below.

classdef intpolys < polygons

  properties
    beta;
  end
  
  methods
    function P = intpolys(varargin)
      P = P@polygons(varargin{:});
      
      P.beta = P.vl;
      for j = 1:P.m
        P.beta(1:P.vc(j),j) = 1 - P.p{j}.alpha;
      end
      P.beta(:,1) = -P.beta(:,1);
    end
    
    function plot(P)
      toHold = ishold;
      hold on
      
      [m vc vl] = polydat(P);
%       xlim = [ inf -inf ]; ylim = xlim;
      for j = 1:m
        w = [ vl(1:vc(j),j); vl(1,j) ];
        plot(real(w),imag(w),'k')
        plot(real(vl(1,j)),imag(vl(1,j)),'ro')
        if j == 1
          xlim = [ min(real(w)) max(real(w)) ];
          ylim = [ min(imag(w)) max(imag(w)) ];
        end
      end
      
      set(gca,'dataaspectratio',[1 1 1]);
      pad = max(0.1*[xlim ylim]);
      axis([xlim(1)-pad xlim(2)+pad ylim(1)-pad ylim(2)+pad])
      
      if ~toHold
        hold off
      end
    end % plot
    
  end %methods
  
end
