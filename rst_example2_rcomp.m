%% resistor example
% irregular poly -> circle -> rectangle
clear
clear fungenbd4_5


NN = 3:10;


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


rb = zeros(1,length(NN)); rt = rb;
for n = 1:length(NN)
  N = NN(n);
  
  fprintf('\nCalculating N = %d\n',N);
  
  % circle domain
  Xc = solve_parambd(P,N,fixed,X0);
  [c r t] = circDomain(Xc,vc);
  
  % rectangle domain
  [PR g gp tr] = rst_solve_rect(c,r,t,vl4,N);
%   zr = preVertices(c,r,tr,PR.vc);
  
  rb(n) = real(PR.vl(3));
  rt(n) = real(PR.vl(4));
end

fprintf('\n');
for n = 1:length(NN)
  fprintf('\t%d \t &   %.8f \t &   %.8f   \\\\\n',NN(n),rb(n),rt(n));
end
fprintf('\n');
