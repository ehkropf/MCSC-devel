%% rst_example2 error data graph
clear

load rst_ex2_data

%%
dXN = diff(XN,1,2);
nXN = max(abs(dXN),[],1)';
dRb = abs(diff(Rb));
dRt = abs(diff(Rt));

%%
figure(1),clf
semilogy(repmat(NN(2:end)',1,4),...
  [rsumN(2:end) dRb dRt nXN],'.-')
set(gca,'xtick',NN)
xlabel('N')

legend('\Sigma r_{\nu j}','|R_b(N)-R_b(N-1)|','|R_t(N)-R_t(N-1)|',...
  '||X_N-X_{N-1}||_\infty')
