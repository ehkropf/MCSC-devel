%% resistor inner connector asymmetric data
clear

NN = 3:12;
npts = 100;


%%
P = polygons(...
  polygon([ 4+2i -1+4i -4-1i 2-4i ]),...
  polygon([ 0.5+0.5i 0.5+2i -1+2i -1+0.5i ]),...
  polygon([ 1.5-0.5i -0.5i -2.3i 1.5-2.3i ])...
);

fixed = [0.4-0.1i; 1.5+0.5i];
m = P.m;
a = 0.4;

X0 = [
   0.281237457555731; 0.274393932701057; -0.018710617633105;
   0.317246953928912; -0.009931383803965; -0.526215401491527;
   1.483990640401263; 3.123926436809910; 4.798533041441599;
   5.134228262769486; 6.686412849512207; 8.230103195309768;
   9.866260107587824; 0.427367537128364; 1.998734689946237;
   3.823310413214275; 5.121931321326818;
];



%%
Rmean = zeros(length(NN),1);
Rmin = Rmean; Rmax = Rmin;
rsumN = Rmax; ctime = Rmax;
XN = zeros(size(X0,1),length(NN));

for n = 1:length(NN)
  N = NN(n);
  
  fprintf('\n\nCalculating for N = %d:\n',N)
  
  tic
  Xc = solve_parambd(P,N,fixed,X0,0);
  ctime(n) = toc;
  XN(:,n) = Xc;
  
  [f C z w rsum] = make_mapbd(Xc,P,N,fixed);
  rsumN(n) = rsum(end);
  
  [c r t] = circDomain(Xc,P.vc);
  zz0 = repmat(c(2:3),npts,1) + ...
    repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,2)*diag(r(2:3));
  [g gp] = map_radslitring(c,r,a,N);
  eta = g(zz0);
  
  rout = abs(eta(:,1)); rin = abs(eta(:,2));
  Rmean(n) = 0.5*log(sum(rout)/sum(rin))/pi;
  Rmin(n) = 0.5*log(min(rout)/max(rin))/pi;
  Rmax(n) = 0.5*log(max(rout)/min(rin))/pi;
end

save rst_ex_inasym_data NN XN rsumN Rmean Rmin Rmax ctime
