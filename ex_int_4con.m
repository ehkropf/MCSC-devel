%% bounded 3 polygon example
clear


%% Polygon container
c = 0.4; d = 0.1; a = 0.2i; b = -0.3;
P = intpolys(...
  [1i; -1.2; -0.8-0.7i; -1i; 1.2+0.3i],...
  [c-d+a; c-d*1i+a+0.07; c+d+a; c+d*1i+a-0.07],...
  [-c+d+b; -c+d*1i+b; -c-d+b; -c-d*1i+b],...
  [ 0.15-0.2i; -0.37-0.34i; -0.17-0.52i ] ...
);

% initial guess
Cg = circdomain(...
  { 0 1 [ 0; pi/2; 2*pi/3; pi; 3*pi/2 ] },...
  { -0.4i 0.2 [ pi/2; pi; 3*pi/2; 0 ] },...
  { 0.4i 0.2 [ 3*pi/2; 0; pi/2; pi ] },...
  { -0.4 0.2 [ 0; 2*pi/3; 4*pi/3 ] } ...
);

fixed = [0.0; 0.0];


%%
opts = intmapopts;
opts.monitor = true;
opts.fignum = 1;

C = circdomain(Cg);
f = intmap(P, fixed, C, opts);



% plot(f)
