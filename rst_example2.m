%% Resistor example.
clear

N = 4;                % reflection level


%%
% L-shape with two rectangular holes.
P = intpolys(...
  polygon([ 3+3i -3+3i -3 0 -3i 3-3i ]), ...
  polygon([ 1 1-1.7i 2-1.7i 2 ]), ... % bottom
  polygon([ 1i 2i -1.5+2i -1.5+1i ]) ...  % top
);
% figure(1),clf,plot(P),return

% Determines f(fixed_pt(1)) = fixed_pt(2).
fixed_pt = [0; 0.7+0.7i];

% The 4 vertices of P that will define the generalized quadralateral.
genquad = [2, 3, 5, 6];


%%
% Circle domain inital guess. From previous trial and error.

% Cg = circdomain(...
%   {0, 1, (0:5)/3*pi}, ...
%   {-0.5i, 0.1, [0.5, 1, 1.5, 2]*pi}, ...
%   {0.5i, 0.1, [1.5, 2, 2.5, 3]*pi} ...
% );

Cg = circdomain(...
  { 0 1 [ 0 0.4766 0.5439 0.9999 1.4565 1.5239 ]*pi },...
  { 0.0160-0.6556i 0.2304 [ 0.6166 1.3437 1.5596 2.0773 ]*pi },...
  { 0.0154+0.6555i 0.2303 [ 1.3835 1.9222 2.4408 2.6566 ]*pi } ...
);


%%
% Find circle domain/get SC map for P.

opts = intmapopts;
opts.monitor = true;
opts.fignum = 1;
opts.N = N;

f = intmap(P, fixed_pt, circdomain(Cg), opts);
C = f.C;

return


%%
% Find rectangle domain given C.

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
