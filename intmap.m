%INTMAP   MCSC interior (bounded) map.
%  This is the interface class for a numerical implementation of the
%  bounded case of the multiply connected Schwarz-Christoffel map. It is a
%  subclass of mcscmap.
%
%  f = intmap(P, F, C) creates a bounded map object, P is a polygons
%  container object, F is a 2-vector of a point in the polygonal and a
%  corresponding point in the circle domain, and C is a circdomain
%  container object which represents the inital guess for the numerical
%  solver.
%
%  f = intmap(P, F, C, opts) uses the intmapopts object to control the
%  behavior of the intmap construction. See intmapopts help for more
%  information.
%
%  intmap methods:
%    map_eval   evaluates the map at given points
%    plot       plot the map
%
%  See also intmapopts, intmap.
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

classdef intmap < mcscmap

  properties
    fixp;
    app_acc;
  end
  
  methods
    function M = intmap(P,fixed,varargin)
      % need to solve the parameter problem here
      % INPUTS:
      %   P = polygons object (polygonal domain)
      %   N = truncation parameter
      %   fixed = interior fixed point: [z,w]
      %   C = circle domain
      %   name-value pairs for extra stuff?
      
      if nargin == 0
        sargs = {};
        fixed = [];
      else
        sargs{1} = P;
        
        Cig = [];
        opts = [];
        for j = 1:numel(varargin)
          if isobject(varargin{j})
            switch class(varargin{j})
              case 'circdomain'
                Cig = varargin{j};
              case 'intmapopts'
                opts = varargin{j};
              otherwise
                error('mcsc:badarg',...
                  'Unrecognized class %s in argument #%d.',...
                  class(varargin{j}),j+2);
            end
          else
            error('mcsc:badarg','Argument #%d unrecognized.',j+2)
          end
        end
        
        if isempty(opts)
          opts = intmapopts;
        end
        if isempty(Cig)
          error('mcsc:noinitguess','Must give a circle domain for now.')
        end
        
        sargs{2} = opts.N;
        sargs{3} = Cig;

        % KLUDGE
        if isempty(opts.objfun)
          sargs{4} = intobjfun(P,fixed,sargs{3},opts);
        else
          sargs{4} = opts.objfun(P,fixed,sargs{3},opts);
        end
        
        if ~opts.nosolver
          sargs{5} = numcontin(sargs{4},sargs{3}.Xu); % solver
          sargs{5}.tolerance = opts.tolerance;
        else
          sargs{5} = [];
        end
      end  
      
      % super class calls the solver and saves the
      % integrand and quadrature objects (each handle classes)
      M = M@mcscmap(sargs{:});
      
      if nargin
        M.fixp = fixed;
        
        M.A = calc_constant(M);
        M.app_acc = calc_apparent_accuracy(M);
        fprintf('Apparent (vertex) accuracy is %.5g.\n\n',M.app_acc);
      end
    end % constructor
    
    function acc = calc_apparent_accuracy(M)
      qobj = M.Q; A = M.A;
      [m vc vl] = polydat(M.P);
      tkj = M.C.t;
      zkj = prevertices(M.C);
      
      w = vl;
      
      % expanded radii
      rr = M.C.r;
      rr(1) = M.C.r(1) - closest_circle(M.C,1)/2;
      for j = 2:m
        rr(j) = M.C.r(j) + closest_circle(M.C,j)/2;
      end
      
      % C_1
      % inward line
      left = zkj(2:vc(1)-1,1);
      right = rr(1)*exp(1i*tkj(2:vc(1)-1,1));
      lpvn = (2:vc(1)-1)'; rpvn = 0*lpvn;
      wo = lineq(qobj, left,lpvn, right,rpvn);
      
      % arc
      left = tkj(2:vc(1)-1,1);
      right = tkj(3:vc(1),1);
      lpvn = 0*left; rpvn = lpvn;
      cv = 0*left; rv = repmat(rr(1),size(left));
      wa = arcq(qobj, left,lpvn, right,rpvn, cv,rv);
      
      % outward line
      left = rr(1)*exp(1i*tkj(3:vc(1),1));
      right = zkj(3:vc(1),1);
      rpvn = (3:vc(1))'; lpvn = 0*rpvn;
      wi = lineq(qobj, left,lpvn, right,rpvn);
      
      w(3:vc(1),1) = vl(2:vc(1)-1,1) + A*(wo + wa + wi);
            
      % other circles
      c = M.C.c;
      for j = 2:m
        idx = sum(vc(1:j-1));
        left = zkj(1,j);
        lpvn = idx + 1;
        w(1,j) = M.fixp(2) - A*lineq(qobj, left,lpvn, M.fixp(1),0);
        
        % outward line
        left = zkj(1:vc(j)-1,j);
        right = c(j) + rr(j)*exp(1i*tkj(1:vc(j)-1,j));
        lpvn = idx + (1:vc(j)-1)';
        rpvn = 0*right;
        wo = lineq(qobj, left,lpvn, right,rpvn);
        
        % arc
        left = tkj(1:vc(j)-1,j);
        right = tkj(2:vc(j),j);
        lpvn = 0*left; rpvn = lpvn;
        cv = repmat(c(j),size(left)); rv = repmat(rr(j),size(left));
        wa = arcq(qobj, left,lpvn, right,rpvn, cv,rv);
        
        % inward line
        left = c(j) + rr(j)*exp(1i*tkj(2:vc(j),j));
        right = zkj(2:vc(j),j);
        rpvn = idx + (2:vc(j))';
        lpvn = 0*rpvn;
        wi = lineq(qobj, left,lpvn, right,rpvn);
        
        w(2:vc(j),j) = vl(1:vc(j)-1,j) + A*(wo + wa + wi);
      end
      
      acc = max(abs(vl(:)-w(:)));
    end % apparent accuracy
    
    function [arc,line,k,j] = integpath(M,z)
      % based on something Driscoll wrote in 2005
      
      vc = M.P.vc; tkj = M.C.t;
      c = M.C.c; r = M.C.r;
      
      [j,t] = nearestcircle(M,z);
      
      % distance (and direction to nearest prevertex)
      dt = tkj(1:vc(j),j) - t;
      dt = mod(dt+pi,2*pi)-pi; % dt in [-pi,pi)
      [mindt,k] = min(abs(dt));
      
      if mindt < 100*eps
        % no arc needed (z is the prevertex or at the same angle on a
        % line to the circle center)
        arc = [];
        line.lpvn = sum(vc(1:j-1)) + k;
      else
        arc.left = tkj(k,j);
        arc.lpvn = sum(vc(1:j-1)) + k;
        arc.right = tkj(k,j) - dt(k);
        line.lpvn = 0;
      end
      
      if (j > 1 && (abs(z-c(j)) < r(j) + 100*eps(r(j)))) ...
          || (j == 1 && abs(z-c(j)) > r(j) + 100*eps(r(j)))
        % call it on the circle
        line = [];
      else
        line.left = c(j)+r(j)*exp(1i*t);
      end
    end % integpath

    function plot(M,opts) %#ok<INUSD>
      if nargin < 2
        % plot options
        pfig = 4;
        cfig = 5;
        sfig = 6;
        npts = 100;
        nrad = 14;
        ncirc = 8;
      else
        error('mcsc:plotopts','Plot options need to be implemented!');
      end
      
      gray = 0.7*ones(1,3);
      
      m = M.P.m;
      c = M.C.c; c = c(:).';
      r = M.C.r; r = r(:)';
      z0 = M.fixp(1);
      
      %% polygons
      figure(pfig)
      clf,hold on
      paxh = gca;
      plot(M.P)
      drawnow
      
      %% circle boundary and slitmap image
      S = intradslitmap(z0,c,r);
      z = repmat(c,npts,1) + ...
        repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,m)*diag(r);
      zet = S(z);
      
      cfh = figure(cfig);
      clf,hold on
      caxh = gca;
      plot(M.C)
%       axis equal
%       axis([-1.1 1.1 -1.1 1.1])
      drawnow
      
      %% polar grid in slit domain
      zetr = repmat(linspace(0,1,npts)',1,nrad) .* ...
        repmat(exp(2i*pi*(0:nrad-1)/nrad),npts,1);
      clrad = linspace(1/(ncirc+1),1-1/(ncirc+1),ncirc);
      zetc = repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,ncirc)*diag(clrad);
            
      figure(sfig)
      clf,hold on
      plot(zetr, 'color', gray)
      plot(zetc, 'color', gray)
      plot(zet,'k')
      plot(real(S(1)),imag(S(1)),'ro')
      axis equal
      axis([-1.1 1.1 -1.1 1.1])
      drawnow
      
      
      %% slit domain geometry info
      tpv = pretip_angles(S);
      vertmod = abs( S(repmat(c(2:end),2,1) + ...
        repmat(r(2:end),2,1).*exp(1i*tpv)) );
      saLR = zeros(size(tpv));
      for j = 2:m
        % slit is probably not exactly at one angle; record deviation.
        % going counterclockwise, call the min angle the left angle, and
        % the max angle the right angle.
        mna = min(mod(angle(zet(:,j)),2*pi));
        mxa = max(mod(angle(zet(:,j)),2*pi));
        saLR(:,j-1) = [ mna; mxa ];
      end
      
      
      %% invert radial lines
      fprintf('\n\nInverting grid lines from radial slit domain...')
      
      zr = zeros(size(zetr));

      % we know where the origin goes...
      zr(1,:) = z0;
      
      % initial guess
      zn = zetr(2,:) + z0;
      
      for k = 2:size(zr,1)
        zn = mapinv_newt(S,zetr(k,:),zn);
        
        plot(caxh,zn,'b.'),drawnow
        
        zr(k,:) = zn;
      end
      
      % redraw circle domain
      set(0,'currentfigure',cfh)
      clf,hold on
      caxh = gca;
      plot(M.C)
      plot(zr, 'color', gray)
%       axis equal
%       axis([-1.1 1.1 -1.1 1.1])
      drawnow
      
      
      %% invert circular lines
      
      % gather info about jumping across radial slits
      % if a point is "too far" from one side of the slit, add one closer
      zg = zeros(size(zetc));
      for j = 2:m
        % find out how many points before the slit (to the left)
        pta = mod(angle(zetc(1:end-1,1)),2*pi);
        b = pta < saLR(1,j-1);
        b = sum(b);
        
        % get the logical indexes in the width of the slit
        sL = vertmod(1,j-1) <= abs(zetc(1,:)) & ...
          abs(zetc(1,:)) <= vertmod(2,j-1);
        
        % add points if needed
        aspc = 0.001; % an aribrary spacing away from circle bdry
        if sum(sL)
          if saLR(1,j-1) - pta(b,1) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(1,j-1)-aspc)); ...
              zetc(b+1:end,:) ];
            zg = [ zg(1:b,:); zg(end,:); zg(b+1:end,:) ];
            b = b + 1;
          end
          if b+1 <= numel(pta) && (pta(b+1) - saLR(2,j-1)) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(2,j-1)+aspc)); ...
              zetc(b+1:end,:) ];
            zg = [ zg(1:b,:); zeros(1,ncirc); zg(b+1:end,:) ];
          elseif pta(1)-2*pi-saLR(2,j-1) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(2,j-1)+aspc)) ];
            zg = [ zg(1:b,:); zeros(1,ncirc) ];
          end
        end
        
        % markers for jumping across slits
        zg(b+1,sL) = j;
      end
      
      zc = zeros(size(zetc));
      
      % initial guess: start with the image of the first radial line...
      % (hope there's no slit here...)
      zetrdx = zeros(size(clrad));
      for k = 1:ncirc
        zetrdx(k) = ...
          find(abs(zetr(:,1)-clrad(k)) == min(abs(zetr(:,1)-clrad(k))),1);
      end
      zn = zr(zetrdx,1).';
      zn = mapinv_odeig(S,zetc(1,:),zn);
      
      for k = 1:size(zc,1)
        wk = zetc(k,:);
        
        % handle slit jumps
        if any(zg(k,:))
          zgi = unique(zg(k,zg(k,:)~=0));
          
          for p = zgi
            zgL = zg(k,:) == p;
            tLE = tpv(1,p-1); tRE = tpv(2,p-1);
            if tLE > tRE, tRE = tRE + 2*pi; end
            
            % a continued line should have a type of symmetry on the circle; we
            % assume the "outgoing" line has about the same angular distance from
            % the slit left-end prevertex as the "incoming" line, up to a scale
            % factor based on the two arc-lengths between the slit prevertices.
            s = 2*pi/(tRE-tLE)-1;
            tzb = mod(angle(zn(zgL) - c(p)),2*pi);
            tzb = tzb + 2*pi*(tLE > tzb) - tLE;
            tzt = tLE - s*tzb;
            % the 1.05 here is a fudge factor
            zn(zgL) = c(p) + 1.05*r(p)*exp(1i*tzt);
          end
        end
        
        zn = mapinv_newt(S,wk,zn);
        
        plot(caxh,zn,'b.'),drawnow
        
        zc(k,:) = zn;
      end
      
      % make cell array for plotting
      zcc = cell(1,sum(sum(zg~=0)+1));
      ci = 1;
      for k = 1:size(zc,2)
        zgi = find(zg(:,k)~=0);
        b = 1;
        for p = 1:numel(zgi)
          zcc{ci} = zc(b:zgi(p)-1,k);
          b = zgi(p);
          ci = ci + 1;
        end
        zcc{ci} = zc(b:end,k);
        ci = ci + 1;
      end
      
      % redraw circle domain
      set(0,'currentfigure',cfh)
      clf,hold on
      plot(zr, 'color', gray)
      for k = 1:numel(zcc)
        plot(zcc{k}, 'color', gray)
      end
      plot(M.C)
%       axis equal
%       axis([-1.1 1.1 -1.1 1.1])
      drawnow
      
      fprintf('done.\n')
      
      
      %% polygonal domain
      fprintf('Evaluating the MCSC map on the grid, be patient...\n')
      
      wr = zr;
      secs = 0;
      fprintf('...evaluating radial grid images...\n')
      for k = 1:nrad
        tic
        wr(:,k) = M.map_eval(zr(:,k));
        secs = secs + toc;
        plot(paxh,wr(:,k), 'color', gray),drawnow
      end
      wc = cell(size(zcc));
      fprintf('...evaluating circular grid images...\n')
      for k = 1:numel(zcc)
        tic
        wc{k} = M.map_eval(zcc{k});
        secs = secs + toc;
        plot(paxh,wc{k}, 'color', gray),drawnow
      end
      fprintf('Map evaluation took %-3.5f seconds.\n\n',secs);
      
    end % plot
    
  end % methods
  
end
