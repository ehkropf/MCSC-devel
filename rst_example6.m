%% resistor example
% irregular poly -> circle -> rectangle
clear

N = 6;
show_prmprb = 1;


%% irregular polygon
P = intpolys(...
  polygon([ 0.5i -0.5+0.5i -0.5+0.75i -1.5+0.75i -1.5-0.5i 0.5-0.5i ...
            0.5-0.75i 1.5-0.75i 1.5+0.5i ]),...
  polygon([ -0.5 -0.75+0.2i -1 -0.75-0.2i ]),...
  polygon([ 0.5 0.75-0.2i 1 0.75+0.2i ])...
);

return
vl4 = [ 3 4 7 8 ];
fixed = [0; 0];
[m vc vl beta] = polydat(P);

% X0 = [...
%   0.2; 0.2;
%   0; 0.5; 0; -0.5;
%   2*pi*(1:8)'/9;
%   2*pi*(0:3)'/4 + 3*pi/2;
%   2*pi*(0:3)'/4 + pi/2;
% ];

X0 = [
  0.1; 0.1;
  0.01; 0.79; -0.01; -0.79;
  1.21; 1.32; 1.45; 1.58; 4.35; 4.45; 4.6; 4.72;
  4.68; 6.63; 7.72; 8.9;
  1.54; 3.49; 4.58; 5.76;
];


%*************************************
%% get the circle domain and polygon map

if show_prmprb && show_prmprb
  figure(1),clf
end

Xc = solve_parambd(P,N,fixed,X0,show_prmprb);

% return

[c r t] = circDomain(Xc,vc);

f = make_mapbd(Xc,P,N,fixed);



%************************************
%% find rectangle domain

[PR g gp tr] = rst_solve_rect(c,r,t,vl4,N);
zr = preVertices(c,r,tr,PR.vc);

% fprintf('\nrN%d = [ %.15f   %.15f ];\n\n',N,real(PR.vl(3:4,1))')


% bottom, then top
% rN7  = [ 2.842004764216455   2.842004764121497 ];
% rN8  = [ 2.841997747152968   2.841997747152820 ];
% rN9  = [ 2.841998946309831   2.841998946309174 ];
% rN10 = [ 2.841998423908737   2.841998423908773 ];
% rN11 = [ 2.841998498461617   2.841998498461600 ];


% return

%************************************
%% plots

figure(2),clf,hold on
plot_rst_poly(P,vl4)

figure(3),clf,hold on
plot_rst_circle(c,r,t,zr,vl4)

figure(4),clf,hold on
plot_rst_rect(PR)

drawnow



%*************************************
%% flowlines
inarc = t(vl4(1:2),1);
opts = rst_flowopts;
opts.npts = 200;
fprintf('\nInverting horizontal lines in rectangle domain ...\n')
tic
[wfr zf] = rst_flowlines(PR,c,zr,inarc,g,gp,opts);
s = toc;
fprintf('\t... took %.4f seconds\n',s);

%% equipotential
barc = t(vl4(2:3),1);
fprintf('\nInverting equipotential lines...\n');
tic
[wpr zpc] = rst_eplines(PR,c,r,tr,barc,g,gp,opts);
s = toc;
fprintf('\t...took %.4f seconds.\n\n',s);

%% polygonal domain
fprintf('Evaluating lines in polygonal domain ...\n');
tic
wf = f(zf);
wpc = cell(size(zpc));
for k = 1:length(zpc)
  wpc{k} = f(zpc{k});
end
s = toc;
fprintf('\t... took %.4f seconds.\n\n',s);



%*******************************
%%
figure(2),clf,hold on
plot(wf,'b')
for k = 1:length(wpc)
  plot(wpc{k},'b')
end
plot_rst_poly(P,vl4)

figure(3),clf,hold on
plot(zf,'b')
for k = 1:length(zpc)
  plot(zpc{k},'b')
end
plot_rst_circle(c,r,t,zr,vl4)

figure(4),clf,hold on
plot(wfr,'b')
plot(wpr,'b')
plot_rst_rect(PR)
