function [wf zf] = rst_flowlines(PR,c,zr,inarc,f,fp,opts)
% generate flow lines in given rectangle domain and return these points and
% their pre-images in the circle domain

% E. Kropf, 2010

if nargin < 7 || isempty(opts)
  opts = rst_flowopts;
end

nlines = opts.nlines;
npts = opts.npts;
newtol = 1e-12;
maxiter = 20;


x = linspace(0,min(real(PR.vl(3:4,1))),npts);
y = linspace(1,0,nlines+2); y = y(2:end-1); 
[x,y] = meshgrid(x,y);
wf = complex(x,y);
zg = wf; zg(:) = nan;


% figure out where grid line are too close to slits; set these to NaN.
% this will have problems if only one point is in the length of a slit
% but I'm lazy right now.
for j = 2:PR.m
  LE = min(real(PR.vl(1:2,j)));
  RE = max(real(PR.vl(1:2,j)));
  CLE = zr(real(PR.vl(1:2,j)) == LE,j);
  CRE = zr(real(PR.vl(1:2,j)) == RE,j);
  msk = abs(imag(wf) - imag(PR.vl(1,j))) < 0.01;
  msk = msk & ( LE <= real(wf) & real(wf) <= RE );
  if nnz(msk)
    mski = find(msk);
    wfi = imag(wf(mski(1)));
    wf(msk) = nan;
    wf(mski([1 length(mski)])) = complex([LE RE],wfi);
    % map hints for keeping preimage points away from circle boundaries
    zg(mski([1 length(mski)])) = ([CLE CRE] - c(j))*1.05 + c(j);
  end
end


%% Newton iteration for flowlines in circle domain
% newtol = 1e-12; maxiter = 20;
zf = zeros(size(wf));

% need an initial guess
td = linspace(inarc(1),inarc(2),nlines+2)';
td = td(2:end-1);
zn1 = exp(1i*td);

for k = 1:size(wf,2)
  if k == 1, zn = zn1; end
  
  wk = wf(:,k); zgk = zg(:,k);
%   done = false(size(zn));
  nanL = isnan(wk);
  zgL = ~isnan(zgk);
  zn(zgL) = zgk(zgL);
%   done(nanL) = true;
  zn(nanL) = nan;
  
%   ni = 0;
%   while ~all(done) && ni < maxiter
%     F = f(zn(~done)) - wk(~done);
%     zn(~done) = zn(~done) - F./fp(zn(~done));
%     done(~done) = abs(F) < newtol;
%     ni = ni + 1;
%   end

  zn = mapinv_newt(wk,f,fp,zn,newtol,maxiter);
  
  zf(:,k) = zn;
end

wf = wf.';
zf = zf.';
