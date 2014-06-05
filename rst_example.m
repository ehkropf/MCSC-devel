%% Circle domain to slit rectangle.
clear
clear fungenbd_rst2


N = 4;
npts = 30;


%%
% the only thing that matters for P here is the length of
% the first side and the turning angles.
P = polygons(...
  polygon([1i 0 3 3+1i]),...
  polygon([2 1]),...
  polygon([2 1]),...
  polygon([2 1])...
);

X0 = [pi; 0; pi; 0; pi; 0];

[m vc vl beta] = polydat(P);

c = [0 -0.25+0.4i 0.4 -0.3-0.4i];
r = [1 0.1 0.2 0.1];
t = zeros(4,m);
t(:,1) = [ 9*pi/8; 11*pi/8; 2*pi; 19*pi/8; ];
t(1:2,2:end) = reshape(X0,2,m-1);
% fixed = [0; 0.75+.45i];



%%**********************************
X0u = unconstr_rst(X0);
% figure(2), clf
tic
[Xu,lam] = contparam(@(X) fungenbd_rst2(X,P,N,c,r,t), X0u);
tm = toc;

fprintf('\nlambda = %1.15f\n',lam);

fprintf('\ninf norm of objective function = %1.15g\n',...
  norm(fungenbd_rst2(Xu,P,N,c,r,t),inf));

fprintf('\nContinuation took %-3.5f seconds.\n\n',tm);

Xc = constr_rst(Xu);
t(1:2,2:end) = reshape(Xc,2,m-1);



% reflect interior circles to exterior
z = preVertices(c,r,t,vc);
cnu = c; rnu = r;
znu = reshape(z,max(vc),1,m);
for j=2:m
  [cnu(1,j),rnu(1,j)] = centrad1(c(1),r(1),c(j),r(j));
  znu(1:vc(j),1,j) = reflect1(c(1),r(1),z(1:vc(j),j));
end
% do N levels of reflection
[znu,cnu,rnu,snu,jlr,rsum] = reflectzsmi(znu,cnu,rnu,cnu,vc,N);


% find induced polygon vertices
w = P.vl;
C = findC(P,cnu,t,znu,N);
ngj = 32;
fp = @(z) fp3bd(z,znu,cnu,beta,vc,N);

w(3,1) = C*gjint(fp, t(2,1),t(3,1), beta(2,1),beta(3,1),...
  ngj,'arc',0,1) + w(2,1);
w(4,1) = C*gjint(fp, t(3,1),t(4,1), beta(3,1),beta(4,1),...
  ngj,'arc',0,1) + w(3,1);

for j = 2:m
  w(1,j) = C*gjint(fp, znu(1,1,1),znu(1,2,j), beta(1,1),beta(1,j),...
    ngj,'line') + w(1,1);
  w(2,j) = C*gjint(fp, t(1,j),t(2,j), beta(1,j),beta(2,j),...
    ngj,'arc',c(j),r(j)) + w(1,j);
end

p = cell(1,m);
for j = 1:m
  p{j} = polygon(w(1:vc(j),j));
end
PR = polygons(p);




%%********************************
% circle domain
figure(1),clf,hold on
eit = exp(2i*pi*(0:100)'/100);
for j = 1:m
  plot(c(j) + r(j)*eit,'k')
end
plot(exp(1i*(t(1,1) + (0:50)*(t(2,1)-t(1,1))/50)),...
  'k','LineWidth',2.5)
plot(exp(1i*(t(3,1) + (0:50)*(t(4,1)-t(3,1))/50)),...
  'k','LineWidth',2.5)
plot(znu(:,1,1),'r*')
for j = 2:m
  plot(znu(1,2,j),'ro')
  plot(znu(2,2,j),'r*')
end
axis equal
axis([-1.2 1.2 -1.2 1.2]);


% slit rectangle domain
figure(2),clf,hold on
plot(w(1:2,1),'k','LineWidth',2.5)
plot(real(w(2:3,1)),imag(w(2:3,1)),'k')
plot(w(3:4,1),'k','LineWidth',2.5)
plot(w([4 1],1),'k')
for j = 2:m
  plot(w(1:2,j),'k')
end
axis equal
ax = axis;
axis([ax(1)-0.1 ax(2)+0.1 ax(3)-0.1 ax(4)+0.1]);

drawnow;



%***********************************
% [Xp,Yp] = meshgrid(linspace(-1,1,npts));
% Zp = complex(Xp,Yp);
% L = abs(Zp) < 1;
% for j = 2:m
%   L = L & (abs(Zp - c(j)) > r(j));
% end
% Zp(~L) = NaN;
% 
% f = mcscmapper(PR,C,N,znu,cnu,rnu,t);
% 
% fprintf('\n\nEvaluating map at grid points...\n');
% tic
% Wp = f(Zp);
% sec = toc;
% fprintf('...took %.4f secs.\n\n',sec);
% 
% figure(1)
% contour(real(Zp),imag(Zp),imag(Wp),9,'b')
% % contour(real(Zp),imag(Zp),real(Wp),9,'b')
% 
% figure(2)
% contour(real(Wp),imag(Wp),imag(Wp),9,'b')
% % contour(real(Wp),imag(Wp),real(Wp),9,'b')

