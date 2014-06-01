%% resistor example
% irregular poly -> circle -> rectangle
clear
clear fungenbd4_5

N = 5;
show_prmprb = 1;


%% irregular polygon
P = polygons(...
  polygon([ -2.8+0.8i -0.9-2i 2.5-1.1i 1.8+1.6i -0.6+2.3i ]),...
  polygon([ -1+0.3i -1.4 -0.9-0.5i -0.4 ]),...
  polygon([ 0.8+0.3i 0.9-0.3i 1.4-0.3i 1.3+0.3i ]),...
  polygon([ 0.1+0.9i 0.7+0.9i 0.2+1.5i ])...
);

vl4 = [ 1 2 3 5 ];
fixed = [0; -0.4+0.6i];
[m vc vl beta] = polydat(P);

% X0 = [...
%   0.2; 0.2; 0.2;
%   0; 0.5; -0.5; 0; 0; -0.5;
%   2*pi*(1:4)'/5;
%   2*pi*(0:3)'/4 + 3*pi/2;
%   2*pi*(0:3)'/4;
%   2*pi*(0:2)'/3 + pi/2;
% ];

X0 = [
  0.2; 0.1; 0.15;
  0.27; 0.33; -0.55; 0.39; -0.45; -0.15;
  1.46; 2.41; 3.06; 4.59;
  4.87; 6.23; 7.96; 9.55;
  5.23; 7.01; 8.38; 9.91;
  0.36; 2.17; 4.54
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

fprintf('\nrN%d = [ %.15f   %.15f ];\n\n',N,real(PR.vl(3:4,1))')

% bottom, then top
rN5 = [ 0.994030176779630   0.993980366165269 ];
rN6 = [ 0.993980109838058   0.993988286395400 ];
rN7 = [ 0.993979802867434   0.993975222162694 ];


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
fprintf('\nInverting horizontal lines in rectangle domain ...\n')
tic
[wfr zf] = rst_flowlines(PR,c,zr,inarc,g,gp);
s = toc;
fprintf('\t... took %.4f seconds\n',s);

%% equipotential
barc = t(vl4(2:3),1);
fprintf('\nInverting equipotential lines...\n');
tic
[wpr zpc] = rst_eplines(PR,c,r,tr,barc,g,gp);
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
