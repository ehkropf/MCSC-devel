%CIRCDOMAIN   MCSC circle domain.
%
%   This handle class is used to represent both bounded and unbounded
%   circle domains. It takes care of constrained/unconstrained
%   transformations as well as sending notifications to listening objects
%   on a change in the domain configuration (see, e.g., one of the
%   objective funcitons).
%
%   C = circdomain(...
%     { 0 1 [ t(1,1); t(2,1); ... ] },...
%     { c_2 r_2 [ t(1,2); t(2,2); ...] },...
%     { c_3 r_3 [ t(1,3); t(2,3); ...] },...
%     ...
%   ) creates a CIRCDOMAIN object representing the circle domain with
%   centers c_1 = 0, c_2, c_3, ...; radii r_1 = 1, r_2, r_3, ...; and
%   prevertex angles t(k,j)
%
%   methods:
%     list them
%
%   See also extobjfun, intobjfun.
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

classdef circdomain < handle

  properties
%     circ; % array of circles
    c;
    r;
    t;
  end
      
  properties (SetObservable)
    Xu;
  end
  
  properties (SetAccess=protected)
    m;
    vc;
  end
 
  properties (Dependent,SetAccess=private)
    Xc;
  end
  
  properties (Access=private)
    internal;
    eit;
  end
  
  events
    externXuUpdate;
  end
  
  methods
    function C = circdomain(varargin)
      if nargin > 0
        C.internal = true;
        
        if isa(varargin{1},'polygon')
          error('mcsc:circdomain','Auto circle domain not implemented.');
        elseif isa(varargin{1},'circdomain')
          O = varargin{1};
          C.m = O.m;
          C.c = O.c;
          C.r = O.r;
          C.vc = O.vc;
          C.t = O.t;
          C.Xu = O.Xu;
        else
          m = numel(varargin);
          C.m = m;
          C.c = zeros(m,1);
          C.r = C.c;
          C.vc = C.c;
          
          t = cell(1,m);
          maxt = 0;
          for j = 1:m
            [C.c(j) C.r(j) t{j}] = deal(varargin{j}{:});
            C.vc(j) = numel(t{j});
            maxt = max(maxt,C.vc(j));
          end
          C.t = zeros(maxt,m);
          for j = 1:m
            C.t(1:C.vc(j),j) = t{j}(:);
          end
        
          C.Xu = unconstrained(C);
        end
        
        C.eit = exp(2i*pi*(0:200)'/200);
        
        C.internal = false;
      end
    end % constructor

    function set.Xu(C,Xu)
      if C.internal
        C.Xu = Xu;
        return
      end
      
      if isequal(C.Xu,Xu)
        return
      end
      
      % unconstrained to constrained transform
      m = C.m; vc = C.vc;
      C.internal = true;

      % radii
      C.r = [1; exp(Xu(1:m-1))];
      
      % copy centers
      cri = Xu(m:3*m-3);
      C.c = [0; complex(cri(1:2:end-1),cri(2:2:end))];
      
      % prevertices (first circle)
      b = 3*m-2;
      sk = cumsum([1;exp(Xu(b+(0:vc(1)-2)))]); %#ok<PROP>
      C.t(2:vc(1),1) = 2*pi*sk(1:end-1)/sk(end); %#ok<PROP>
%       phi1 = 2*pi/(1 + sum(exp(Xu(b:b+vc(1)-2)))); %#ok<PROP>
%       for k = 2:vc(1) %#ok<PROP>
%         C.t(k,1) = phi1*(1+sum(exp(Xu(b:b+k-3))));
%       end
      
      % prevertices on other circles
      b = b + vc(1)-1; %#ok<PROP>
      for j = 2:m
        C.t(1,j) = mod(Xu(b),2*pi);
        sk = cumsum([1;exp(Xu(b+(1:vc(j)-1)))]); %#ok<PROP>
        C.t(2:vc(j),j) = C.t(1,j) + 2*pi*sk(1:end-1)/sk(end); %#ok<PROP>
%         phi1 = 2*pi/(1 + sum(exp(Xu(b+1:b+vc(j)-1)))); %#ok<PROP>
%         for k = 2:vc(j) %#ok<PROP>
%           C.t(k,j) = C.t(1,j) + phi1*(1+sum(exp(Xu(b+1:b+k-2))));
%         end
        b = b + vc(j); %#ok<PROP>
      end
      C.internal = false;
      
      C.Xu = Xu;
      
      notify(C,'externXuUpdate');
    end

    function set.c(C,c)
      if C.internal
        C.c = c;
        return
      end
      
      if numel(c) ~= C.m
        error('MCSC:wrongsize','Number of elements must match connectivity.')
      end
      
      C.c = c(:);
      C.internal = true;
      C.Xu = unconstrained(C);
      C.internal = false;
    end
    
    function set.r(C,r)
      if C.internal
        C.r = r;
        return
      end
      
      if numel(r) ~= C.m
        error('MCSC:wronsize','Number of elements must match connectivity.')
      end
      
      C.r = r(:);
      C.internal = true;
      C.Xu = unconstrained(C);
      C.internal = false;
    end
    
    function set.t(C,t)
      if C.internal
        C.t = t;
        return
      end
      
      if size(t) ~= size(C.t)
        error('MCSC:wrongsize','Check theta matrix; wrong size.')
      end
      
      for j = 1:C.m
        C.t(1:C.vc(j),j) = t(1:C.vc(j),j);
      end
      C.internal = true;
      C.Xu = unconstrained(C);
      C.internal = false;
    end
    
    function Xc = get.Xc(C)
      Xc = constrained(C);
    end
    
    function Xu = unconstrained(C)
      % constrained to unconstrained transform
      m = C.m; vc = C.vc;
      Xu = zeros(3*m-4+sum(C.vc),1);
      
      % radii
      Xu(1:m-1) = log(C.r(2:m)); %log(Xc(1:m-1));
      
      % centers don't change
      Xu(m:2:3*m-4) = real(C.c(2:m));
      Xu(m+1:2:3*m-3) = imag(C.c(2:m));
      
      % prevertices (first circle)
      b = 3*m-2;
      phi = diff([0; C.t(2:vc(1),1); 2*pi]); %#ok<PROP>
      phi = mod(phi,2*pi); % enforce sum(phi)=2*pi (needed?)
      Xu(b+(0:vc(1)-2)) = log(phi(2:vc(1))/phi(1)); %#ok<PROP>
      
      % prevertices on other circles
      b = b + vc(1)-1; %#ok<PROP>
      for j = 2:m
        Xu(b) = C.t(1,j);
        phi = diff([C.t(1:vc(j),j); Xu(b)+2*pi]); %#ok<PROP>
        phi = mod(phi,2*pi);
        Xu(b+(1:vc(j)-1)) = log(phi(2:vc(j))/phi(1)); %#ok<PROP>
        b = b + vc(j); %#ok<PROP>
      end
      
      if sum(isinf(Xu))
        warning('mcsc:unconstr',['At least one of the unconstrained '...
          'variables is at infinity.'])
      end
      
%       C.Xu = Xu;
    end
    
    function Xc = constrained(C)
      m = C.m; vc = C.vc;
      
      Xc = zeros(sum(C.vc)+3*m-4,1);
      Xc(1:m-1) = C.r(2:end);
      Xc(m-1+(1:2:2*(m-1)-1)) = real(C.c(2:end));
      Xc(m-1+(2:2:2*(m-1))) = imag(C.c(2:end));
      
      b = 3*m-2;
      Xc(b+(0:vc(1)-2)) = C.t(2:vc(1),1); %#ok<PROP>
      b = b + vc(1)-1; %#ok<PROP>
      for j = 2:m
        Xc(b+(0:vc(j)-1)) = C.t(1:vc(j),j); %#ok<PROP>
        b = b + vc(j); %#ok<PROP>
      end
    end
        
    function zkj = prevertices(C)
      zkj = C.t;
      
      for j = 1:C.m
        zkj(1:C.vc(j),j) = C.c(j) + C.r(j)*exp(1i*C.t(1:C.vc(j),j));
      end
    end % prevertices
    
    function plot(C)
%       ah = gca;
      toHold = ishold;
      
      cbox = [inf -inf inf -inf];
      zkj = prevertices(C);
      hold on
      
      for j = 1:C.m
        plot(C.c(j) + C.r(j)*C.eit,'k-')
        cbox([1 3]) = ...
          min(cbox([1 3]),[real(C.c(j)) imag(C.c(j))]-C.r(j));
        cbox([2 4]) = ...
          max(cbox([2 4]),[real(C.c(j)) imag(C.c(j))]+C.r(j));
        plot(real(zkj(1,j)),imag(zkj(1,j)),'ro')
        plot(real(zkj(2:C.vc(j),j)),imag(zkj(2:C.vc(j),j)),'r*')
      end
      cbox(1:2) = mean(cbox(1:2)) + 0.6*diff(cbox(1:2))*[-1 1];
      cbox(3:4) = mean(cbox(3:4)) + 0.6*diff(cbox(3:4))*[-1 1];
      set(gca,'dataaspectratio',[1 1 1])
      axis(cbox)
%       axis([-1.1 1.1 -1.1 1.1])
      
      if ~toHold, hold off, end
    end
    
    function replicate(C,n)
      m = C.m; c = C.c; r = C.r; t = C.t/pi;
      if nargin < 2 || isempty(n)
        n = 4;
      end
      
      fprintf('\nC = circdomain(...\n')
      for j = 1:m
        if j ~= 1
          cr = sprintf('  { %%0.%df%%+0.%dfi %%0.%df [',n,n,n);
          cr = sprintf(cr,real(c(j)),imag(c(j)),r(j)); %#ok<PROP>
          nt = C.vc(j);
          vn = 1:nt;
        else
          cr = sprintf('  { 0 1 [ 0');
          nt = C.vc(j)-1;
          vn = 2:nt+1;
        end
        pverts = repmat(' %%0.%df',1,nt);
        pverts = sprintf(pverts,repmat(n,1,nt));
        pverts = sprintf(pverts,t(vn,j)); %#ok<PROP>
        if j ~= m
          dots = ' ]*pi },...\n';
        else
          dots = ' ]*pi } ...\n';
        end
        
        fprintf([cr pverts dots])
      end
      fprintf(');\n\n')
    end % replicate
    
    function [d,p] = closest_circle(C,j)
      % give the index p of the circle closest to C_j and the distance d
      % between them
      
      if all(abs(C.c(2:end)) < C.r(1))
        % internal case
        if j == 1
          [d,p] = min( C.r(1) - abs(C.c(2:end)) - C.r(2:end) );
        else
          [d,p] = min( [ C.r(1) - abs(C.c(j)) - C.r(j); ...
            abs(C.c([false (2:C.m~=j)]) - C.c(j)) ...
              - C.r([false (2:C.m~=j)]) - C.r(j) ] );
        end
        
      else
        % external case
        [d,p] = min( abs(C.c(1:C.m~=j) - C.c(j)) ...
          - C.r(1:C.m~=j) - C.r(j) );
      end
      
      if p >= j
        p = p + 1;
      end
    end % closest_circle
    
  end % methods
    
end
