%% resistor interior contacts, nonsymmetric
clear
clear fungenbd4_5

N = 6;
colorcirc = 0;


%%
% P = polygons(...
%   polygon([ 4 4+4i -4+4i -4-4i 4-4i ]),...
%   polygon([ -0.5+0.5i -0.5+2i -2+2i -2+0.5i ]),...
%   polygon([ 1.5-0.5i -0.5i -2.3i 1.5-2.3i ])...
% );

P = polygons(...
  polygon([ 4+2i -1+4i -4-1i 2-4i ]),...
  polygon([ 0.5+0.5i 0.5+2i -1+2i -1+0.5i ]),...
  polygon([ 1.5-0.5i -0.5i -2.3i 1.5-2.3i ])...
);

% fixed = [0.3+0.2i; 1+1i];
fixed = [0.4-0.1i; 1.5+0.5i];
m = P.m;

% X0 = [
%   0.2; 0.2;
%   -0.2; 0.3; 0.2; -0.3;
%   pi/2; pi; 3*pi/2;
%   7*pi/4; pi/4; 3*pi/4; 5*pi/4;
%   pi/4; 3*pi/4; 5*pi/4; 7*pi/4;
% ];


X0 = [
   0.281237457555731; 0.274393932701057; -0.018710617633105;
   0.317246953928912; -0.009931383803965; -0.526215401491527;
   1.483990640401263; 3.123926436809910; 4.798533041441599;
   5.134228262769486; 6.686412849512207; 8.230103195309768;
   9.866260107587824; 0.427367537128364; 1.998734689946237;
   3.823310413214275; 5.121931321326818;
];
% Xc = X0;


%%
% figure(1),clf
Xc = solve_parambd(P,N,fixed,X0,1);

% return

[c r t] = circDomain(Xc,P.vc);
f = make_mapbd(Xc,P,N,fixed);

figure(2),clf,hold on

for j = 1:m
  verts = P.vl([1:P.vc(j) 1],j);
  if j == 1
    plot(verts,'k')
  else
    plot(verts,'k','linewidth',2.5)
  end
end
set(gca,'dataaspect',[1 1 1])



a = 0.4;

npts = 100;
zz0 = repmat(c,npts,1) + ...
  repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,m)*diag(r);

figure(3),clf,hold on
plot(zz0)
if colorcirc && colorcirc
  plot(zz0)
else
  for j = 1:m
    if j == 1
      plot(zz0(:,j),'k')
    else
      plot(zz0(:,j),'k','linewidth',2.5)
    end
  end
end
axis equal


[g gp] = map_radslitring(c,r,a,N);
eta = g(zz0);

figure(4),clf,hold on
if colorcirc && colorcirc
  plot(eta)
else
  for j = 1:m
    if j == 1
      plot(eta(:,j),'k')
    else
      plot(eta(:,j),'k','linewidth',2.5)
    end
  end
end
axis equal

resistance = 0.5*log(sum(abs(eta(:,2)))/sum(abs(eta(:,3))))/pi;
fprintf('\n  resistance = %.15f\n\n',resistance)

pause
% return

%% flowlines in slitring domain
rlines = 13;
rmax = min(abs(eta(:,2)));
rmin = max(abs(eta(:,3)));

npts = 100;
rl = logspace(log10(rmin),log10(rmax),npts)';
rl = repmat(rl,1,rlines)*diag(exp(2i*pi*(0:rlines-1)/rlines));

figure(4)
if colorcirc && colorcirc
  plot(rl,'k')
else
  plot(rl,'b')
end

drawnow




%% invert flowlines
zn = c(3) + r(3)*exp(2i*pi*(0:rlines-1)'/rlines);
wa = angle(g(zn(1)));
zn = c(3) + (zn-c(3))*exp(1i*-wa);


rl = rl.';
z = zeros(size(rl));
zn = mapinv_odeig(rl(:,1),g,gp,zn);

for k = 1:size(rl,2)
  zn = mapinv_newt(rl(:,k),g,gp,zn);
  z(:,k) = zn;
  figure(3),plot(zn,'kx'),drawnow
end
z = z.';

figure(3),clf,hold on
if colorcirc && colorcirc
  plot(zz0)
  plot(z,'k')
else
  for j = 1:m
    if j == 1
      plot(zz0(:,j),'k')
    else
      plot(zz0(:,j),'k','linewidth',2.5)
    end
  end
  plot(z,'b')
end
axis equal

drawnow



%% equipotential lines in slitring domain
eplines = 8;
oa0 = max(angle(eta(:,1)));
oa1 = min(angle(eta(:,1)));
% crad = linspace(rmin,rmax,eplines+2);
crad = logspace(log10(rmin),log10(rmax),eplines+2);
crad = crad(2:end-1);

cl = exp(1i*((2*pi-(abs(oa0-oa1)))*(0:npts-1)'/(npts-1) + oa0));
cl = repmat(cl,1,eplines)*diag(crad);

figure(4)
if colorcirc && colorcirc
  plot(cl,'k')
else
  plot(cl,'b')
end

drawnow


%% invert circles
cl = cl.';

zn = exp(1i*(7*pi/8*(0:eplines-1)/(eplines-1) - pi/2));
zn = mapinv_odeig(cl(:,1),g,gp,zn);

zc = zeros(size(cl));
for k = 1:size(cl,2)
  zn = mapinv_newt(cl(:,k),g,gp,zn);
  zc(:,k) = zn;
  figure(3),plot(zn,'kx'),drawnow
end
zc = zc.';

figure(3),clf,hold on
if colorcirc && colorcirc
  plot(zz0)
  plot(z,'k')
  plot(zc,'k')
else
  for j = 1:m
    if j == 1
      plot(zz0(:,j),'k')
    else
      plot(zz0(:,j),'k','linewidth',2.5)
    end
  end
  plot(z,'b')
  plot(zc,'b')
end
axis equal

drawnow



%%
fprintf('\nEvaluating lines in polygonal domain...\n')
tic
w = f(z);
wc = f(zc);
s = toc;
fprintf('\t...took %f seconds.\n\n',s)

figure(2)
plot(w,'b')
plot(wc,'b')
