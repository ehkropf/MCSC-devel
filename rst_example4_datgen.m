%% resistor example
% irregular poly -> circle -> rectangle
clear
clear fungenbd4_5

% N = 15;
show_prmprb = 0;
NN = 3:6;


%% irregular polygon
P = polygons(...
  polygon([ 0 2i -2+2i -2-2i 2-2i 2 ]),...
  polygon([ -0.75 -0.08 ])...
);
vl4 = [ 2 3 4 5 ];
fixed = [-0.3i; -0.75-0.75i];
[m vc vl beta] = polydat(P);

X0 = [...
  0.17;
  0.47; 0.18;
  1.32; 1.5; 4.06; 5.28; 5.36;
  2.94; 6.05
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
