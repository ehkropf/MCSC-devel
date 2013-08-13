%EXTMAP   MCSC exterior (unbounded) map.
%  This is the interface class for a numerical implementation of the
%  unbounded case of the multiply connected Schwarz-Christoffel map. It is
%  a subclass of mcscmap.
%
%  f = extmap(P, C) creates an unbounded map object, P is a polygons
%  container object and C is a circdomain container object which represents
%  the initial guess for the numerical solver.
%
%  f = extmap(P, C, opts) uses the extmapopts object to control the
%  behavior of the extmap construction. See extmapopts help for more
%  information.
%
%  extmap methods:
%    map_eval   evalutes the map at given points
%    plot       plot the map
%
%  See also extmapopts, mcscmap.
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

classdef extmap < mcscmap

  properties
    apparent_accuracy;
    app_acc;
  end
  
  methods
    function M = extmap(P,varargin)
      if nargin == 0
        sargs = {};
      else
        sargs{1} = P;
        
        Cig = [];
        opts = [];
        for j = 1:numel(varargin)
          if isobject(varargin{j})
            switch class(varargin{j})
              case 'circdomain'
                Cig = varargin{j};
              case 'extmapopts'
                opts = varargin{j};
              otherwise
                error('mcsc:badarg',...
                  'Unrecognized class %s.',class(varargin{j}));
            end
          else
            error('mcsc:badarg','Unrecognized argument.')
          end
        end
        
        if isempty(opts)
          opts = extmapopts;
        end
        if isempty(Cig)
          error('mcsc:noinitguess','Must give a circle domain for now.')
        end
        
        sargs{2} = opts.N;
        sargs{3} = Cig;
                        
        sargs{4} = extobjfun(P,sargs{3},opts); % objective function
        
        if ~opts.nosolver
          % FIXME: solver construction here should be generalized!
          % (too tightly bound to numcontin right now)
          sargs{5} = numcontin(sargs{4},sargs{3}.Xu); % solver
          sargs{5}.hmin = opts.numcontin.hmin;
          
          sargs{5}.tolerance = opts.tolerance;
        else
          sargs{5} = [];
        end
      end
      
      % super class calls the solver and saves the
      % integrand and quadrature objects (each handle classes)
      M = M@mcscmap(sargs{:});
      
      if nargin
        M.A = calc_constant(M);
        M.apparent_accuracy = calc_apparent_accuracy(M);
        fprintf('Apparent (vertex) accuracy is %.5g.\n', ...
          M.apparent_accuracy);
      end
      
    end % constructor
    
    function acc = calc_apparent_accuracy(M)
      qobj = M.Q; A = M.A;
      c = M.C.c; r = M.C.r; tkj = M.C.t;
      zkj = prevertices(M.C);
      [m vc vl] = polydat(M.P);
      
      w = vl;
      
      % expanded radii
      rr = r;
      for j = 1:m
        rr(j) = r(j) + closest_circle(M.C,j)/2;
      end
            
      % C_1
      % outward line
      left = zkj(2:vc(1)-1,1);
      right = rr(1)*exp(1i*tkj(2:vc(1)-1,1));
      lpvn = (2:vc(1)-1)'; rpvn = 0*lpvn;
      wo = lineq(qobj, left,lpvn, right,rpvn);
        
      % arc
      left = tkj(2:vc(1)-1,1);
      right = tkj(3:vc(1),1);
      lpvn = 0*left; rpvn = lpvn;
      cv = 0*left; rv = rr(1)*ones(size(left));
      wa = arcq(qobj, left,lpvn, right,rpvn, cv,rv);
      
      % inward line
      left = rr(1)*exp(1i*tkj(3:vc(1),1));
      right = zkj(3:vc(1),1);
      rpvn = (3:vc(1))'; lpvn = 0*rpvn;
      wi = lineq(qobj, left,lpvn, right,rpvn);
      
      w(3:vc(1),1) = vl(2:vc(1)-1,1) + A*(wo + wa + wi);
      
      
      % other circles
      for j = 2:m
        idx = sum(vc(1:j-1));
        lpvn = idx + 1; left = zkj(1,j);
        w(1,j) = w(1,1) - A*lineq(qobj, left,lpvn, zkj(1,1),1);
        
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
        cv = c(j)*ones(size(left)); rv = rr(j)*ones(size(left));
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
    
    function val = calc_sidelen_compare(M)
      qobj = M.Q; A = M.A;
      [m vc vl] = polydat(M.P);
      tkj = M.C.t; c = M.C.c; r = M.C.r;
      
      w = vl;
      
      for j = 1:m
        off = sum(vc(1:j-1));
        lpvn = off + (1:vc(j))'; rpvn = [lpvn(2:end); off+1];
        left = tkj(1:vc(j),j); right = [left(2:end); tkj(1,j)];
        right(left>right) = right(left>right) + 2*pi;
        cv = repmat(c(j),vc(j),1); rv = repmat(r(j),vc(j),1);
        
        w(1:vc(j),j) = ...
          A*arcq(qobj, left,lpvn, right,rpvn, cv,rv);
      end
      
      vld = zeros(sum(vc),1);
      wd = vld;
      for j = 1:m
        b = sum(vc(1:j-1));
        vld(b+(1:vc(j))) = abs(diff([vl(1:vc(j),j); vl(1,j)]));
        wd(b+(1:vc(j))) = abs(w(1:vc(j),j));
      end
      
      val = max(abs(wd(:)-vld(:)));      
    end % side-length comparison
    
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
      
      if abs(z-c(j)) < r(j) + 100*eps(r(j))
        % call it on the circle
        line = [];
      else
        line.left = c(j)+r(j)*exp(1i*t);
      end
    end % integpath

    function plot(M, z0, opts) %#ok<INUSD>
      if nargin < 3
        % plot options
        pfig = 4;
        cfig = 5;
        sfig = 6;
        npts = 200;
        nrad = 15;
        ncirc = 11;
      else
        error('mcsc:plotopts','Plot options need to be implemented!');
      end
      
      gray = [0 0 1]; %0.7*ones(1,3);
      
      m = M.P.m;
      c = M.C.c; c = c(:).';
      r = M.C.r; r = r(:)';
      
      %% polygons
      figure(pfig)
      clf,hold on
      paxh = gca;
      plot(M.P)
      drawnow
      
      %% circle boundary and slitmap image
      z = repmat(c,npts,1) + ...
        repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,m)*diag(r);
      
%       cbox = [inf -inf inf -inf];
%       for j = 1:m
%         cbox([1 3]) = min(cbox([1 3]),[real(c(j)) imag(c(j))]-r(j));
%         cbox([2 4]) = max(cbox([2 4]),[real(c(j)) imag(c(j))]+r(j));
%       end
%       cbox(1:2) = mean(cbox(1:2)) + diff(cbox(1:2))*[-0.7 0.7];
%       cbox(3:4) = mean(cbox(3:4)) + diff(cbox(3:4))*[-0.7 0.7];
      
      cfh = figure(cfig);
      clf,hold on
      caxh = gca;
      plot(M.C)
%       axis equal
%       axis(cbox)
      drawnow

      if nargin < 2 || isempty(z0)
        fprintf('\nSelect conformal center...');
        [x y] = ginput(1);
        z0 = complex(x,y);
        fprintf('\n\tchosen at x=%.4g y=%.4g.\n\n',x,y);
      end
      
      S = extradslitmap(z0,c,r);
      zet = S(z);

      %% polar grid in slit domain
      maxrad = 1.1*max(abs(zet(:)));
      zetr = repmat(linspace(0,maxrad,npts)',1,nrad) .* ...
        repmat(exp(2i*pi*(0:nrad-1)/nrad),npts,1);
      clrad = maxrad*linspace(1/(ncirc+1),1,ncirc);
      zetc = repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,ncirc)*diag(clrad);
            
      figure(sfig)
      clf,hold on
      plot(zetr,'color', gray)
      plot(zetc,'color', gray)
      plot(zet,'k')
%       plot(real(S(1)),imag(S(1)),'ro')
      axis equal
      axis(maxrad*[-1.1 1.1 -1.1 1.1])
      drawnow
      
      
      %% slit domain geometry information
      tpv = pretip_angles(S);
      vertmod = abs( S(repmat(c,2,1) + repmat(r,2,1).*exp(1i*tpv)) );
      saLR = zeros(size(tpv)); % slit angle deviation
      for j = 1:m
        % slit is probably not exactly at one angle; record deviation.
        % going counterclockwise, call the min angle the left angle, and
        % the max angle the right angle.
        mna = min(mod(angle(zet(:,j)),2*pi));
        mxa = max(mod(angle(zet(:,j)),2*pi));
        saLR(:,j) = [ mna; mxa ];
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
        
        plot(caxh,zn,'.','color',gray),drawnow
        
        zr(k,:) = zn;
      end
      
      % redraw circle domain
      set(0,'currentfigure',cfh)
      clf,hold on
      caxh = gca;
      plot(M.C)
      plot(zr,'color', gray)
      drawnow
      
      
      %% invert circular lines
      
      % gather info about jumping across radial slits
      % if a point is "too far" from one side of the slit, add one closer
      zg = zeros(size(zetc));
      for j = 1:m
        % find out how many points before the slit (to the left)
        pta = mod(angle(zetc(1:end-1,1)),2*pi);
        b = pta < saLR(1,j);
        b = sum(b);
        
        % get the logical indexes in the width of the slit
        sL = vertmod(1,j) <= abs(zetc(1,:)) & ...
          abs(zetc(1,:)) <= vertmod(2,j);
        
        % add points if needed
        aspc = 0.001; % an aribrary spacing away from circle bdry
        if sum(sL)
          if saLR(1,j) - pta(b,1) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(1,j)-aspc)); ...
              zetc(b+1:end,:) ];
            zg = [ zg(1:b,:); zeros(1,ncirc); zg(b+1:end,:) ];
            b = b + 1;
          end
          if b+1 <= numel(pta) && (pta(b+1) - saLR(2,j)) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(2,j)+aspc)); ...
              zetc(b+1:end,:) ];
            zg = [ zg(1:b,:); zg(end,:); zg(b+1:end,:) ];
          elseif pta(1)-2*pi-saLR(2,j) > aspc
            zetc = [ zetc(1:b,:); zetc(1,:)*exp(1i*(saLR(2,j)+aspc)) ];
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
          %     zgL = zg(k,:) ~= 0;
          zgi = unique(zg(k,zg(k,:)~=0));
          
          for p = zgi
            zgL = zg(k,:) == p;
            tLE = tpv(1,p); tRE = tpv(2,p);
            if tLE > tRE, tRE = tRE + 2*pi; end
            
            % a continued line should have a type of symmetry on the circle; we
            % assume the "outgoing" line has about the same angular distance from
            % the slit left-end prevertex as the "incoming" line, up to a scale
            % factor based on the two arc-lengths between the slit prevertices.
            s = 2*pi/(tRE-tLE)-1;
            tzb = mod(angle(zn(zgL) - c(p)),2*pi);
            tzb = tzb + 2*pi*(tLE > tzb) - tLE;
            tzt = tLE - s*tzb;
            % the 1.05 here is a kludge factor (stay away from actual
            % boundary!)
            zn(zgL) = c(p) + 1.05*r(p)*exp(1i*tzt);
          end
        end
        
        zn = mapinv_newt(S,wk,zn);
        
        plot(caxh,zn,'.','color',gray),drawnow
        
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
      plot(zr,'color',gray)
      for k = 1:numel(zcc)
        plot(zcc{k},'color',gray)
      end
      plot(M.C)
%       axis equal
%       axis(cbox)
      drawnow
      
      fprintf('done.\n')
      
      
      %% polygonal domain
      fprintf('Evaluating the MCSC map on the grid, be patient...\n')
      
      wr = zr;
      secs = 0;
      fprintf('...evaluating radial lines...\n')
      for k = 1:nrad
        tic
        wr(:,k) = M.map_eval(zr(:,k));
        secs = secs + toc;
        plot(paxh,wr(:,k),'color',gray),drawnow
      end
      wc = cell(size(zcc));
      fprintf('...evaluating circular lines...\n')
      for k = 1:numel(zcc)
        tic
        wc{k} = M.map_eval(zcc{k});
        secs = secs + toc;
        plot(paxh,wc{k},'color',gray),drawnow
      end
      fprintf('Map evaluation took %-3.5f seconds.\n\n',secs);
      
    end % plot
    
    function plot_factor(M, k, j)
      c = M.C.c; r = M.C.r; m = M.C.m;
      
      n = 100;
      zc = repmat(c(:).',n,1) ...
        + repmat(exp(2i*pi*(0:n-1)'/(n-1)),1,m)*diag(r);
      
      figure(7),clf,hold on
%       cm = colormap(lines);
%       style = {'-','--','-.',':'};
%       for p = 1:m
%         plot(zc(:,p),'linestyle',style{mod(p-1,4)+1},'color',cm(p,:))
%       end
      plot(zc)
      zkj = c(j) + r(j)*exp(1i*M.C.t(k,j));
      plot(real(zkj),imag(zkj),'ko')
      
      cbox([1 3]) = min([real(c(:)) imag(c(:))]-[r(:) r(:)]);
      cbox([2 4]) = max([real(c(:)) imag(c(:))]+[r(:) r(:)]);
      cbox(1:2) = mean(cbox(1:2)) + 0.55*diff(cbox(1:2))*[-1 1];
      cbox(3:4) = mean(cbox(3:4)) + 0.55*diff(cbox(3:4))*[-1 1];
      set(gca,'dataaspectratio',[1 1 1])
      axis(cbox)
      
      figure(8),clf,hold on
      wc = M.fp(zc,k,j);
%       for p = 1:m
%         plot(wc(:,p),'linestyle',style{mod(p-1,4)+1},'color',cm(p,:))
%       end
      plot(wc)
      wkj = M.fp(zkj,k,j);
      plot(real(wkj),imag(wkj),'ko')
      plot(1,0,'kx')
      set(gca,'dataaspectratio',[1 1 1])
      axis([-0.1 2.1 -1.1 1.1])
    end % plot_factor
    
    %%%%% temporary property name migration code
    function set.app_acc(M, value)
      M.apparent_accuracy = value; %#ok<MCSUP>
    end
    
    function aa = get.app_acc(M)
      aa = M.apparent_accuracy;
    end
      
  end % methods
  
end
