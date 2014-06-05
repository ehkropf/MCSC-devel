%% Resistor example.
clear

N = 4;                % reflection level


%%
% L-shape with two rectangular holes.
P = intpolys(...
  polygon([ 3+3i -3+3i -3 0 -3i 3-3i ]), ...
  polygon([ 1.8-0.3i 1.2-0.3i 1.2-1.8i 1.8-1.8i ]), ...
  polygon([ -0.2+1.7i -1.9+1.7i -1.9+1.2i -0.2+1.2i ]) ...
);

% P = intpolys(...
%   polygon([ 3+3i -3+3i -3 0 -3i 3-3i ]), ...
%   polygon([ 1.25 1.25-0.5i 1.75-0.5i 1.75 ]), ...
%   polygon([ 1.25i 1.75i -1+1.75i -1+1.25i ]) ...
% );
figure(1),clf,plot(P),return

% Determines f(fixed_pt(1)) = fixed_pt(2).
fixed_pt = [0; 0.8+0.8i];


%%
% Circle domain inital guess. From previous trial and error.

% Cg = circdomain(...
%   {0, 1, [0, 1.59, 1.78, 3.13, 4.47, 4.68]}, ...
%   {-0.13-0.69i, 0.16, [1.19, 2.53, 4.37, 5.04]}, ...
%   {-0.12+0.71i, 0.15, [5.16, 7.43, 8.18, 9.93]} ...
% );

Cg = circdomain(...
  {0, 1, (0:5)/3*pi}, ...
  {-0.5i, 0.1, [0.5, 1, 1.5, 2]*pi}, ...
  {0.5i, 0.1, [1.5, 2, 2.5, 3]*pi} ...
);

% Cg = circdomain(...
%   {0, 1, (0:5)/3*pi}, ...
%   {0.5i, 0.1, [0.75; 0; 0.25; 0.5]*2*pi}, ...
%   {-0.5i, 0.1, [0.25; 0.5; 0.75; 1]*2*pi} ...
% );

% return
% vl4 = [ 5 6 3 4 ];
% [m vc vl beta] = polydat(P);
% 
% X0 = [...
%   0.16; 0.15;
%   -0.13; -0.69; -0.12; 0.71;
%   1.59; 1.78; 3.13; 4.47; 4.68;
%   1.19; 2.53; 4.37; 5.04;
%   5.16; 7.43; 8.18; 9.93;
% ];


%%
% Find circle domain/get SC map.

opts = intmapopts;
opts.monitor = true;
opts.fignum = 1;
opts.N = N;
f = intmap(P, fixed_pt, circdomain(Cg), opts);

% Xc = solve_parambd(P,N,fixed,X0,1);
% [c r t] = circDomain(Xc,vc);
% 
% f = make_mapbd(Xc,P,N,fixed);

return

%************************************
%% find rectangle domain

[PR g gp tr] = rst_solve_rect(c,r,t,vl4,N);
zr = preVertices(c,r,tr,PR.vc);

% bottom, then top here...
rN9  = [ 1.919740023295317   1.919732662327783 ];
rN10 = [ 1.919740009719658   1.919739482637953 ];
rN11 = [ 1.919740862425454   1.919739959498808 ];
rN12 = [ 1.919740781104729   1.919740709299187 ];

return

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
% flowlines
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



%************************************
%% plots

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
