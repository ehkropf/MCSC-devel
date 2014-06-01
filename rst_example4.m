%% resistor example
% irregular poly -> circle -> rectangle
clear
clear fungenbd4_5

N = 15;
show_prmprb = 1;


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



%*************************************
%% get the circle domain and polygon map

if show_prmprb && show_prmprb
  figure(1),clf
end
Xc = solve_parambd(P,N,fixed,X0,show_prmprb);
[c r t] = circDomain(Xc,vc);

f = make_mapbd(Xc,P,N,fixed);



%************************************
%% find rectangle domain

[PR g gp tr] = rst_solve_rect(c,r,t,vl4,N);
zr = preVertices(c,r,tr,PR.vc);

fprintf('\nrN%d = [ %.15f   %.15f ];\n\n',N,real(PR.vl(3:4,1))')

% bottom, then top here
rN9  = [ 1.843548140064338   1.843548736410529 ];
rN10 = [ 1.843548224567550   1.843548224567531 ];
rN11 = [ 1.843548233414226   1.843548264408984 ];
rN12 = [ 1.843548237806245   1.843548237806266 ];
rN13 = [ 1.843548238266055   1.843548239876998 ];
rN14 = [ 1.843548238494318   1.843548238494305 ];
rN15 = [ 1.843548238518247   1.843548238601994 ];


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
flopts = rst_flowopts;
flopts.npts = 200;
fprintf('\nInverting horizontal lines in rectangle domain ...\n')
tic
[wfr zf] = rst_flowlines(PR,c,zr,inarc,g,gp,flopts);
s = toc;
fprintf('\t... took %.4f seconds\n',s);

%% equipotential
barc = t(vl4(2:3),1);
fprintf('\nInverting equipotential lines...\n');
tic
[wpr zpc] = rst_eplines(PR,c,r,tr,barc,g,gp,flopts);
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
