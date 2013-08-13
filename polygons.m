%POLYGONS   MCSC polygon container superclass.
%   This handle class is the parent of the polygonal domain subclasses.
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

classdef polygons < handle

  properties
    p; % cell array of polygons
    m; % connectivity
    vc; % vertex count
    vl; % vertex list
  end
  
  properties (Abstract)
    beta; % subclass must define; external, bounded, etc.
  end
  
  methods
    function P = polygons(varargin)
      if nargin == 0
        P.p = {};
        P.m = 0;
        P.vc = [];
        P.vl = [];
      elseif isa(varargin{1},'polygons')
        P = varargin{1};
      else
        P.m = numel(varargin);
        P.vc = zeros(1,P.m);
        
        if isa(varargin{1},'polygon')
          pl = varargin;
        else
          pl = cell(1,P.m);
          for j = 1:P.m
            pl{j} = polygon(varargin{j});
          end
        end
        
        P.p = pl;
        P.vc = zeros(1,P.m);
        for j = 1:P.m
          P.vc(j) = numel(pl{j}.vertex);
        end
        P.vl = zeros(max(P.vc),P.m);
        for j = 1:P.m
          P.vl(1:P.vc(j),j) = pl{j}.vertex(:);
        end
      end
    end % constructor
    
%     function b = subsref(a,s)
%       badsyntax = 0;
%       switch s(1).type
%         case '{}'
%           b = a.p{s.subs{:}};
%         case '.'
%           if length(s)>1
%             switch s(2).type
%               case '{}'
%                 b = a.(s(1).subs){s(2).subs{:}};
%               case '()'
%                 b = a.(s(1).subs)(s(2).subs{:});
%               otherwise
%                 badsyntax = 1;
%             end
%           else
%             b = a.(s(1).subs);
%           end
%         otherwise
%           badsyntax = 1; 
%       end
%       if badsyntax
%         error('mcsc:badsyntax','Bad use of polygons object.');
%       end
%       
%     end % subsref
    
    function display(P)
      fprintf('\nPolygons in this container:\n')
      for j = 1:P.m
        display(P.p{j})
      end
    end % disp
    
%     function plot(P)
%       toHold = ishold;
%       
%       plot(P.p{1})
%       hold on
%       for j = 2:P.m
%         plot(P.p{j})
%       end
%       
%       if ~toHold
%         hold off
%       end
%     end % plot
    
    function [m,vc,vl,beta] = polydat(P)
      m = P.m;
      vc = P.vc; vl = P.vl;
      beta = P.beta;
    end
    
  end % methods
  
  methods (Abstract)
    plot(P)
  end % abstract methods
end
