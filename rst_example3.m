%% resistor example
% irregular poly -> circle -> rectangle
clear
clear fungenbd4_5


N = 4;
npts = 100;


%% irregular polygon
P = polygons(...
  polygon([ 3+3i -3+3i -3 0 -3i 3-3i ]),...
  polygon([ -0.2+1.2i -0.2+1.7i -1.9+1.7i -1.9+1.2i ]),...
  polygon([ 1.5 0.5 ])...
);
vl4 = [ 2 3 5 6 ];
fixed = [0; 1+1i];
[m vc vl beta] = polydat(P);

X0 = [...
  0.15; 0.11;
  -0.32; 0.63; -0.31; -0.38;
  1.95; 2.15; 3.26; 4.24; 4.35;
  4.71; 6.09; 8.06; 8.72;
  5.7; 9.01;
];



%*************************************
%% get the circle domain and polygon map

% figure(1),clf
Xc = solve_parambd(P,N,fixed,X0,0);
[c r t] = circDomain(Xc,vc);

f = make_mapbd(Xc,P,N,fixed);



%************************************
%% find rectangle domain

[PR g gp tr] = rst_solve_rect(c,r,t,vl4,N);
zr = preVertices(c,r,tr,PR.vc);


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
flopts = rst_flowopts;
flopts.npts = 100;
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
