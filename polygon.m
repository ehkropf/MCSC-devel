%POLYGON   MCSC polygon class.
%
%   p = polygon(w) creates a POLYGON object where w is a list of ordered
%   vertices such that the domain is on the left of the polygon.
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

classdef polygon

  properties
    vertex;
    alpha;
  end
  
  methods
    function p = polygon(w)
      if nargin == 0
        w = [];
      end
      
      if isempty(w)
        p.vertex = [];
        p.alpha = [];
      elseif isa('w','polygon')
        p = w;
      else
        [p.alpha flip] = p.calc_angles(w(:));
        if flip
          p.vertex = flipud(w(:));
          p.alpha = flipud(p.alpha);
        else
          p.vertex = w(:);
        end
      end
    end % constructor
    
    function d = subsref(p,s)
      switch s(1).type
        case '()'
          d = p.vertex(s.subs{:});
        case '.'
          if length(s)>1
            d = p.(s(1).subs)(s(2).subs{:});
          else
            d = p.(s.subs);
          end
        otherwise
          error('Specify value for x as obj(x)')
      end
    end
    
    function display(p)
      fprintf(' Vertex                                    Alpha\n');
      fprintf(' ---------------------------------------   -----------------\n');
      for j = 1:numel(p.vertex)
        if imag(p.vertex(j)) < 0
          sgn = '-';
        else
          sgn = '+';
        end
        if ~real(p.vertex(j))
          fprintf(' %18d %s %17.15fi   %17.15f\n',...
            0,sgn,abs(imag(p.vertex(j))),p.alpha(j));
        elseif ~imag(p.vertex(j))
          fprintf(' %18.15f %s 0i%16c   %17.15f\n',...
            abs(real(p.vertex(j))),sgn,' ',p.alpha(j));
        else
          fprintf(' %18.15f %s %17.15fi   %17.15f\n',...
            real(p.vertex(j)),sgn,abs(imag(p.vertex(j))),p.alpha(j));
        end
      end
      fprintf('\n');
    end
    
    function n = length(p)
      n = length(p.vertex);
    end
    
    function plot(p)
      plot(p.vertex([1:end 1]))
    end %plot
    
  end % methods
  
  methods (Static)
    function [a,flip] = calc_angles(w)
      % from sctoolbox, eventually should just use that code as it does
      % more...
      n = numel(w);
      in = w - w([n 1:n-1]);
      out = in([2:n 1]);
      a = mod( angle(-in.*conj(out))/pi, 2);

      idx = sum(1-a);
      if abs(idx-round(idx)) > 100*sqrt(numel(w))*eps
        error('mcsc:badpoly','Invalid polygon.')
      end
      if sum(1-a) < 0
        flip = true;
      else
        flip = false;
      end
    end    
    
  end % static methods
  
end
