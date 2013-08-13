%% sample bounded domain driver
clear


%% create polygon domain
P = extpolys(...
  [ 0 -1.5-2i 1.5-2i ],...
  [ -1+1i -1.5+(1i+sqrt(3)*1i/2) -2+1i ],...
  [ 1+1i 2+1i 1.5+1i+sqrt(3)*1i/2 ] ...
);


%% circle domain initial guess
Cg = circdomain(...
  { 0 1 [ 0; 2*pi/3; 4*pi/3 ] },...
  { 2.4+1.3i 0.37 [ 3*pi/2; pi/4; 3*pi/4 ] },...
  { 2.4-1.3i 0.37 [ pi/2; 5*pi/4; 7*pi/4 ] } ...
);


%% map options
opts = extmapopts;

%- reflection method is the default
% opts.method = 'refl';

%- finite product/least squares
opts.method = 'fpls';

%- series truncation level is connectivity*N
% opts.N = 25; % default is 25

%- set tolerance or number of Gauss-Jacobi points
%- where log10(opts.tolerance) = opts.ngj
%- since GJ accuracy (tolerance) ~ ngj
% opts.tolerance = 1e-14;
% opts.ngj = 14;

%- visually monitor solver using given figure number
opts.monitor = true;
% opts.fignum = 1; % default is 1


%% create map
%- circdomain is a handle class; clone the initial guess here using the
%- copy constructor
f = extmap(P, circdomain(Cg), opts);
%- solution circle domain
C = f.C;


%% plotting
%- will prompt for selection of conformal center
% plot(f) 

%- Can also specify desired conformal center
% cc = 2.162+0.06212i;
% plot(f, cc)
