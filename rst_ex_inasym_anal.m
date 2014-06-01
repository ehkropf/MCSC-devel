%% data analysis of asymmeric interior example
clear

load rst_ex_inasym_data

dXN = diff(XN,1,2);
nXN = max(abs(dXN),[],1)';
RN = abs(diff(Rmean));
Nd = NN(2:end);

%%
figure(1),clf
semilogy(repmat(Nd',1,3),[rsumN(2:end) RN nXN],'.-')
set(gca,'xtick',NN)
xlabel('N')
legend('\Sigma r_{\nu j}','|R_{mean}(N) - R_{mean}(N-1)|',...
  '||X_N - X_{N-1}||_\infty')
