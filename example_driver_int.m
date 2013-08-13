%% sample bounded domain driver
clear


%% create polygon domain
c=0.4; d=0.1; a=0.2i; b=-0.3;
P = intpolys(...
  [ 1i; -1.2; -0.8-0.7i; -1i; 1.2+0.3i ],...
  [ c-d+a; c-d*1i+a+0.07; c+d+a; c+d*1i+a-0.07 ],...
  [ -c+d+b; -c+d*1i+b; -c-d+b; -c-d*1i+b ] ...
);

%- fixed = [zeta; omega];
%- specifies f(zeta) = omega
fixed = [0; 0];
           
%% circle domain initial guess
Cg = circdomain(...
  { 0 1 [ 0 pi/2 2*pi/3 pi 3*pi/2 ] },...
  { -0.4i 0.2 [ pi/2 pi 3*pi/2 0 ] },...
  { 0.4i 0.2 [ 3*pi/2 0 pi/2 pi ] } ...
);


%% map options
opts = intmapopts;

%- reflection method is the only method currently available for bounded
%- maps; it's also the default
% opts.method = 'refl';

%- reflection truncation level
% opts.N = 6; % default is 6

%- set tolerance or number of Gauss-Jacobi points
%- where log10(opts.tolerance) = opts.ngj
%- since GJ accuracy (tolerance) ~ ngj
% opts.tolerance = 1e-6;
% opts.ngj = 6;

%- visually monitor solver using given figure number
opts.monitor = true;
% opts.fignum = 1; % default is 1


%% create map
%- circdomain is a handle class; clone the initial guess here using the
%- copy constructor
f = intmap(P, fixed, circdomain(Cg), opts);
%- solution circle domain
C = f.C;


%% plotting
% plot(f)
