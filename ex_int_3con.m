%% test crowding
clear


%% Polygon container
c=0.4; d=0.1; a=0.2i; b=-0.3;
P = intpolys(...
  [1i; -1.2; -0.8-0.7i; -1i; 1.2+0.3i],...
  [c-d+a; c-d*1i+a+0.07; c+d+a; c+d*1i+a-0.07],...
  [-c+d+b; -c+d*1i+b; -c-d+b; -c-d*1i+b] ...
);

% initial guess
Cg = circdomain(...
  { 0 1 [ 0; 1/2; 2/3; 1; 3/2 ]*pi },...
  { -0.4i 0.2 [ 1/2; 1; 3/2; 0 ]*pi },...
  { 0.4i 0.2 [ 3/2; 0; 1/2; 1 ]*pi } ...
);

fixed = [0.0; 0.0];


%%
opts = intmapopts;
opts.monitor = true;
opts.N = 4;

C = circdomain(Cg);
f = intmap(P, fixed, C, opts);


% plot(f)
