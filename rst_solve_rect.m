function [PR g gp tr] = rst_solve_rect(c,r,t,vl4,N)
% solve parameter problem for slit rectangle domain with given circle
% domain. returns the slit rectangle domain object, the map g and its
% derivative gp, and the prevertices in the circle domain.

% E. Kropf, 2010

m = length(c);

pr = cell(1,m);
pr{1} = polygon([ 1i 0 1 1+1i ]);
for j = 2:m
  pr{j} = polygon([1 2]);
end
PR = polygons(pr);

% Get the generalized quadralateral vertices.
tr = zeros(4,m);
tr(1,1) = t(vl4(1),1);
for k = 2:4
  tr(k,1) = t(vl4(k),1);
  if tr(k,1) < tr(k-1,1)
    tr(k,1) = tr(k,1) + 2*pi;
  end
end

% Initial circle domain guess.
XR0 = zeros(2*(m-1),1);
for j = 2:m
  XR0(2*j-3:2*j-2) = [tr(1,1),tr(1,1)+pi];
end

XR0u = unconstr_rst(XR0);

tic
[XRu,lam] = contparam(@(X) fungenbd_rst2(X,PR,N,c,r,tr), XR0u);
tm = toc;
fprintf('\nlambda = %1.15f\n',lam);
fprintf('\ninf norm of objective function = %1.15g\n',...
  norm(fungenbd_rst2(XRu,PR,N,c,r,tr),inf));
fprintf('\nContinuation took %-3.5f seconds.\n\n',tm);

XRc = constr_rst(XRu);
tr(1:2,2:end) = reshape(XRc,2,m-1);

z = preVertices(c,r,tr,PR.vc);
cnu = c; rnu = r;
znu = reshape(z,max(PR.vc),1,m);
for j=2:m
  [cnu(1,j),rnu(1,j)] = centrad1(c(1),r(1),c(j),r(j));
  znu(1:PR.vc(j),1,j) = reflect1(c(1),r(1),z(1:PR.vc(j),j));
end
% do N levels of reflection
[znu,cnu,rnu] = reflectzsmi(znu,cnu,rnu,cnu,PR.vc,N);


% find induced polygon vertices
w = PR.vl;
C = findC(PR,cnu,tr,znu,N);
ngj = 32;
fp = @(z) fp3bd(z,znu,cnu,PR.beta,PR.vc,N);

w(3,1) = C*gjint(fp, tr(2,1),tr(3,1), PR.beta(2,1),PR.beta(3,1),...
  ngj,'arc',0,1) + w(2,1);
w(4,1) = C*gjint(fp, tr(3,1),tr(4,1), PR.beta(3,1),PR.beta(4,1),...
  ngj,'arc',0,1) + w(3,1);

for j = 2:m
  if j == 2
    w(1,j) = C*gjint(fp, z(1,1),z(1,j), PR.beta(1,1),PR.beta(1,j),...
      ngj,'line') + w(1,1);
  else
    w(1,j) = C*gjint(fp, z(2,j-1),z(1,j), PR.beta(2,j-1),PR.beta(1,j),...
      ngj,'line') + w(2,j-1);
  end
  w(2,j) = C*gjint(fp, tr(1,j),tr(2,j), PR.beta(1,j),PR.beta(2,j),...
    ngj,'arc',c(j),r(j)) + w(1,j);
end

for j = 1:m
  pr{j} = polygon(w(1:PR.vc(j),j));
end
PR = polygons(pr);



g = mcscmapper(PR,C,N,znu,cnu,rnu,tr);
gp = @(z) C*fp3bd(z,znu,cnu,PR.beta,PR.vc,N);
