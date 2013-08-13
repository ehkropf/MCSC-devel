%% unbounded 10 polygon example
clear


%%
P = extpolys(...
  polygon([0 -1+1i -1-1i]),...
  polygon([-.17-1.8i -.7-2.2i -.44-2.7i .032-2.3i]),...
  polygon([.21-2.8i -.2-3.8i .53-4.2i .5-3.2i]),...
  polygon([.68-1.2i .5-1.7i .97-2.6i 1.3-2.4i 1.3-1.6i]),...
  polygon([.62-.48i 1.5-.8i 1.5-.19i]),...
  polygon([2.2+.48i 2-.4i 2.9-.015i]),...
  polygon([.5+1.1i 1.1+.75i .62+.42i 1.8+.83i 1.6+1.3i]),...
  polygon([.3+1.7i 1+2i 1+2.6i .59+2.4i]),...
  polygon([.41+3.1i .94+3.5i -.29+3.6i]),...
  polygon([-.52+1.7i -.41+2.8i -1.1+2.5i -1.5+2i -.93+2.2i]) ...
);

Cg = circdomain(...
  { 0 1 [ 0 .6 1.4 ]*pi },...
  { 0.4-3.2i 0.5 [ .4 .9 1.4 1.9 ]*pi },...
  { 1.3-5i 0.7 [ .5 1.1 1.7 2.2 ]*pi },...
  { 2.4-2.8i 0.7 [ .6 1 1.5 1.8 2.2 ]*pi },...
  { 2.8-0.8i 0.5 [ 1 1.7 2.3 ]*pi },...
  { 4.4-0.04i 0.5 [ .6 1.3 2 ]*pi },...
  { 2.7+1.3i 0.7 [ 0.9 1 1.2 1.9 2.3 ]*pi },...
  { 2+3.1i 0.5 [ 1.3 1.8 2.3 2.7 ]*pi },...
  { 1.4+5i 0.6 [ 1.5 2.1 2.9 ]*pi },...
  { -0.2+3.2i 0.7 [ 1.7 2.3 2.7 3.1 3.4 ]*pi } ...
);


%%
opts = extmapopts;
opts.method = 'fpls';
opts.N = 10;
opts.monitor = 1;
opts.fignum = 1;

f = extmap(P, circdomain(Cg), opts);


% %%
% cc = 1.57-1.241i;
% plot(f,cc)
