%% resistor bar example
% see Trefethen, "Analysis and design of polygonal resistors by conformal
% mapping." 1984.

clear

fprime = 'refl';
N = 10;
colorcirc = 0;


%% bar-bell resistor
% p.698: "... all side lengths equal to 1, 1/2, or 1/4."
P = polygons(...
  polygon([ 0.25i -0.5+0.25i -0.5+0.5i -1.5+0.5i -1.5-0.5i ...
           -0.5-0.5i -0.5-0.25i 0.5-0.25i 0.5-0.5i 1.5-0.5i ...
           1.5+0.5i 0.5+0.5i 0.5+0.25i ]),...
  polygon([ -0.75 -0.75+0.25i -1.25+0.25i -1.25-0.25i -0.75-0.25i ]),...
  polygon([ 0.75 0.75-0.25i 1.25-0.25i 1.25+0.25i 0.75+0.25i ])...
);

[m vc vl beta] = polydat(P);
fixed = [0; 0];

% X0 = [
%   0.15; 0.15;
%   0; 0.3; 0; -0.3;
%   2*pi*(1:12)'/13;
%   3*pi/2 + 2*pi*(0:4)'/5;
%   pi/2 + 2*pi*(0:4)'/5;
% ];

X0 = [
   0.015086285340753; 0.015086267236461; -0.000166410519941;
   0.981624373360132; -0.000165795426656; -0.981624394245061;
   1.505892839263334; 1.539921029816123; 1.566343547709987;
   1.575590629565584; 1.602012635240045; 1.636039522502695;
   4.647145975676726; 4.681172812113288; 4.707594858922135;
   4.716841910053967; 4.743264213227991; 4.777292251033749;
   4.712716123001660; 6.400725960045659; 7.580050650064988;
   8.128381052774674; 9.307793339843695; 1.570507598767029;
   3.258592824151452; 4.437980801084158; 4.986310618112451;
   6.165658607077876;
];
% Xc = X0;


%
figure(1),clf
sopts.monitor = 1;
sopts.halfrule = 0;
sopts.fprime = fprime;
Xc = solve_param_b(P,fixed,N,X0,sopts);
% Xc = solve_parambd(P,N,fixed,X0,1);

% return

[c r t] = circDomain(Xc,vc);
% f = make_mapbd(Xc,P,N,fixed);
f = make_map_mcsc_b(Xc,P,fixed,N,fprime,sopts);

% return

figure(2),clf,hold on
for j = 1:m
  verts = P.vl([1:P.vc(j) 1],j);
  if j == 1
    plot(verts,'k')
  else
    plot(verts,'k','linewidth',2.5)
  end
end
axis equal
axis([-1.7 1.7 -0.7 0.7])


%%
a = 0.2i;
[g gp] = map_radslitring(c,r,a,N);

npts = 100;
zz = repmat(c,npts,1) + ...
  repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,m)*diag(r);

figure(3),clf,hold on
if colorcirc && colorcirc
  plot(zz)
else
  for j = 1:m
    if j == 1
      plot(zz(:,j),'k')
    else
      plot(zz(:,j),'k','linewidth',2.5)
    end
  end
end
axis equal

eta = g(zz);
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

% rN8  = 2.767848366224551;
% rN9  = 2.768019239814318;
% rN10 = 2.768998783120371;
% rN11 = 2.769055540165757;
% rN12 = 2.768790779067001;

% return

%% flowlines
rlines = 15;
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


% initial initial guess for inversion
zn = [ c(3) + r(3)*exp(1i*(2*pi/3*(0:7)'/7 + pi/2));
       c(3) + r(3)*exp(1i*(-2*pi/3*(7:-1:1)'/7 + pi/2)) ];
% wa = angle(g(zn(1)));
% zn = c(3) + (zn-c(3))*exp(1i*-wa);


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
  plot(zz)
  plot(z,'k')
else
  for j = 1:m
    if j == 1
      plot(zz(:,j),'k')
    else
      plot(zz(:,j),'k','linewidth',2.5)
    end
  end
  plot(z,'b')
end
axis equal

drawnow



%% equipotential lines
eplines = 10;
oa0 = max(mod(angle(eta(:,1)),2*pi));
oa1 = min(mod(angle(eta(:,1)),2*pi));
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

zn = exp(1i*(pi*(1:eplines)/(eplines+1) - pi/2));
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
  plot(zz)
  plot(z,'k')
  plot(zc,'k')
else
  for j = 1:m
    if j == 1
      plot(zz(:,j),'k')
    else
      plot(zz(:,j),'k','linewidth',2.5)
    end
  end
  plot(z,'b')
  plot(zc,'b')
end
axis equal

drawnow




%%
fprintf('\nMapping lines to polygonal domain...\n')
tic
w = f(z);
wc = f(zc);
s = toc;
fprintf('\t...took %f seconds.\n\n',s);

figure(2)
plot(w,'b')
plot(wc,'b')
