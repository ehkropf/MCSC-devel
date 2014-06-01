function [wp zpc] = rst_eplines(PR,c,r,tr,barc,f,fp,opts)
% equipotential lines for resistor problem

% E. Kropf, 2010

if nargin < 8 || isempty(opts)
  opts = rst_flowopts;
end

nlines = opts.nlines;
npts = opts.npts;
newtol = opts.newtol;
maxiter = opts.maxiter;

m = length(c);

x = linspace(0,min(real(PR.vl(3:4,1))),nlines+2);
x = x(2:end-1);
y = linspace(0,1,npts);
[x,y] = meshgrid(x,y);
wp = complex(x,y).';
zg = zeros(size(wp));


% how to handle jumping across slits? awkwardly?
zk = 1;
for j = 2:m
  sh = imag(PR.vl(1,j));
  b = imag(wp(1,:)) < sh;
  b = sum(b);
  LE = min(real(PR.vl(1:2,j)));
  RE = max(real(PR.vl(1:2,j)));
  tLE = tr(real(PR.vl(1:2,j)) == LE,j);
  tRE = tr(real(PR.vl(1:2,j)) == RE,j);
  sL= LE <= wp(:,1) & wp(:,1) <= RE;

  if sum(sL)
    if abs(imag(wp(1,b)) - sh) > 0.01
      wp = [ wp(:,1:b) wp(:,1)+1i*(sh-0.01) wp(:,b+1:end) ];
      zg = [ zg(:,1:b) zg(:,end) zg(:,b+1:end) ];
      b = b+1;
    end
    if abs(imag(wp(1,b+1)) - sh) > 0.01
      wp = [ wp(:,1:b) wp(:,1)+1i*(sh+0.01) wp(:,b+1:end) ];
      zg = [ zg(:,1:b) zg(:,end) zg(:,b+1:end) ];
    end
    
    % hint for initial guess on other side of slit
    zg(sL,b+1) = zk;
    zgh(:,zk) = [c(j) r(j) tLE tRE]; %#ok<AGROW>
    zk = zk+1;
  end
end




%% Newton iteration to trace lines in circle domain
zp = zeros(size(wp));

% initial guess
if barc(2) < barc(1), barc(1) = barc(1)-2*pi; end
td = linspace(barc(1),barc(2),nlines+2)';
td = td(2:end-1);
zn = exp(1i*td);
zn = mapinv_odeig(wp(:,1),f,fp,zn);


for k = 1:size(wp,2)
  wk = wp(:,k);

  % handle jumps across slits using hints from above
  zgL = zg(:,k) > 0;
  if any(zgL)
    zgi = unique(zg(zgL,k));
    for p = 1:length(zgi)
      zgLp = zg(:,k) == zgi(p);
      cg = zgh(1,zgi(p)); rg = zgh(2,zgi(p));
      tLE = zgh(3,zgi(p)); tRE = zgh(4,zgi(p));
      if tLE > pi, tLE = tLE - 2*pi; end
      tRE = mod(tRE,2*pi);
      if tLE > tRE, tRE = tRE + 2*pi; end
      
      % a continued line should have a type of symmetry on the circle; we
      % assume the "outgoing" line has about the same angular distance from
      % the slit left-end prevertex as the "incoming" line, up to a scale
      % factor based on the two arc-lengths between the slit prevertices.
      s = 2*pi/(tRE-tLE)-1;
      tzb = angle(zn(zgLp) - cg) - tLE;
      tzt = -s*tzb + tLE;
      zn(zgLp) = cg + 1.05*rg*exp(1i*tzt);
    end
  end
    
  zn = mapinv_newt(wk,f,fp,zn,newtol,maxiter);
  
  zp(:,k) = zn;
end


wp = wp.';
zp = zp.';
zg = zg';

ci = 1;
for k = 1:size(zp,2);
  zgi = find(zg(:,k)~=0);
  b = 1;
  for p = 1:length(zgi)
    zpc{ci} = zp(b:zgi(p)-1,k); %#ok<AGROW>
    b = zgi(p);
    ci = ci+1;
  end
  zpc{ci} = zp(b:end,k); %#ok<AGROW>
  ci = ci+1;
end
