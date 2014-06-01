%% rst_example2 data generation for analysis
clear
clear fungenbd4_5

% N = 8;
NN = 3:12;

%% irregular polygon
P = polygons(...
  polygon([ 0 -3i 3-3i 3+3i -3+3i -3 ]),...
  polygon([ -0.2+1.2i -0.2+1.7i -1.9+1.7i -1.9+1.2i ]),...
  polygon([ 1.2-0.3i 1.2-1.8i 1.8-1.8i 1.8-0.3i ])...
);
vl4 = [ 5 6 3 4 ];
fixed = [0; 0.6+0.6i];
[m vc vl beta] = polydat(P);

X0 = [...
  0.16; 0.15;
  -0.13; -0.69; -0.12; 0.71;
  1.59; 1.78; 3.13; 4.47; 4.68;
  1.19; 2.53; 4.37; 5.04;
  5.16; 7.43; 8.18; 9.93;
];


%%
Rb = zeros(numel(NN),1);
Rt = Rb;
% Rmean = zeros(length(NN),1);
% Rmin = Rmean; Rmax = Rmean;
rsumN = zeros(numel(NN),1); % rsumSRo = rsumSC; rsumSRi = rsumSC;
XN = zeros(size(X0,1),numel(NN));

for n = 1:length(NN);
  N = NN(n);
  
  fprintf('\n\nCalculating for N = %d:\n',N)
  
  Xc = solve_parambd(P,N,fixed,X0);
  XN(:,n) = Xc;
  
  [f C z w rsum] = make_mapbd(Xc,P,N,fixed);
  [c r t] = circDomain(Xc,vc);
  
  PR = rst_solve_rect(c,r,t,vl4,N);
  Rb(n) = real(PR.vl(3,1));
  Rt(n) = real(PR.vl(4,1));
  
%   a = 0.2i;
%   [g gp rsumgo rsumgi] = map_radslitring(c,r,a,N);
%   
%   zz = repmat(c(2:3),npts,1) + ...
%     repmat(exp(2i*pi*(0:npts-1)'/(npts-1)),1,2)*diag(r(2:3));
%   eta = g(zz);
% 
%   rout = abs(eta(:,1)); rin = abs(eta(:,2));
%   Rmean(n) = 0.5*log(sum(rout)/sum(rin))/pi;
%   Rmin(n) = 0.5*log(min(rout)/max(rin))/pi;
%   Rmax(n) = 0.5*log(max(rout)/min(rin))/pi;
  
  rsumN(n) = rsum(end);
%   rsumSRo(n) = rsumgo(end);
%   rsumSRi(n) = rsumgi(end);
end

save rst_ex2_data NN Rb Rt XN rsumN
