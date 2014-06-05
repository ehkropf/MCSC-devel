% MCSC abstract reflection method class
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
% Parts of this file are due to Tom DeLillo, Wichita State, 2003.

classdef fprefl < integrand

properties (Abstract)
  C;
end

methods (Static)
  function [co, ro] = centrad(c,r,ci,ri)
    % reflects circle with center ci, radius ri through circle with
    % center c, radius r, to get circle with center co, radius ro
    
    co = c + r^2*(ci - c)/(abs(ci - c)^2 - ri^2);
    ro = r^2*ri/abs(abs(ci - c)^2 - ri^2);
  end
  
  function zo = reflect(c,r,zi)
    % reflects zi through circle with center c, radius r to get zo
    zo = c + r^2./conj(zi - c);
  end
  
  function [zn,cn,rn,sn,jlr,rsum] = reflectzsmi(zin,cn,rn,sn,vc,N)
    % Uses method of images (mi), that is,
    % *** ONLY reflections through original m circles are performed ***
    % This code reflects circles exterior to each other
    % with centers, radii cnu, rnu
    % and center reflections snu and reflections of
    % corner preimages znu
    % through each other N times
    % and draws circles and computes and estimates areas
    %
    %   ***** INPUT  (nu =1) *** OUTPUT (nu => 1) ****
    % N = number of (levels of) reflections of circles and corner preimages
    % cnu(nu,j) = center of reflection nu of circle j
    % snu(nu,j) = reflection nu of center from circle j
    % rnu(nu,j) = radius of reflection nu of circle j
    % znu(k,nu,j) = reflection nu of corner k from circle j
    % mj(j) = number of corner preimages on circle j
    % beta(k,j) = kth turning angle on jth circle
    %
    % *** Auxiliary array:  jlr(nu,j) = leading index of snu(nu,j), etc.
    %         i.e., jlr(nu,j) = index of original circle through which
    %         snu was last reflected - indices cannot be repeated.
    %
    
    % Tom DeLillo, Wichita State U., 6/03
    
    % preallocate for speed!
    m = numel(rn);
    zrows = zeros(sum((m-1).^(1:N)),m);
    cn = [cn; zrows];
    rn = [rn; zrows];
    sn = [sn; zrows];
    zn = zeros(max(vc),sum((m-1).^(0:N)),m);
    zn(:,1,:) = zin;
    sumrl = zeros(N,m);
    
    % initialize (nu = 1 is level 0 - original circles)
    % jlr(nu,j) = leading index = index of circle of last reflection
    
%     for j=1:m
%       jlr(1,j)=j; %#ok<AGROW>
%     end
    jlr = [1:m; zrows];
    
    
    %  outer j loop over m circles
    for j=1:m
      %  track reflections from circle j down to level N
      num = 0;
      for level = 1:N
        sumrl(level,j)=0;
        nul = num+1;
        if m ~= 2
          num = ((m-1)^level - 1)/(m-2);
        elseif m == 2
          num = nul;
        end
        nuj = num;
        for nu = nul:num
          for j1 = 1:m
            if j1 ~= jlr(nu,j)
              nuj = nuj+1;
              jlr(nuj,j) = j1;
              [cn(nuj,j), rn(nuj,j)] = ...
                fpintrefl.centrad(cn(1,j1),rn(1,j1),cn(nu,j),rn(nu,j));
              zn(1:vc(j),nuj,j) = ...
                fpintrefl.reflect(cn(1,j1),rn(1,j1),zn(1:vc(j),nu,j));
              sn(nuj,j) = ...
                fpintrefl.reflect(cn(1,j1),rn(1,j1),sn(nu,j));
              sumrl(level,j) = sumrl(level,j) + rn(nuj,j);
            end
          end
        end
      end
    end
    
    rsum = zeros(1,N+1);
    rsum(1) = sum(rn(1,:));
    for level = 1:N
      rsum(level+1) = 0;
      for j = 1:m
        rsum(level+1) = rsum(level+1) + sumrl(level,j);
      end
    end
  end % reflectzsmi
    
end

end
