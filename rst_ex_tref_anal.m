%% data analysis of Trefethen's example
clear

load rst_ex_tref_data

dXN = diff(XN,1,2);
nXN = [ nan; max(abs(dXN),[],1)' ];

Rt = 2.76886750270;
Rmed = (Rmax+Rmin)/2;


%%
figure(1),clf
semilogy(repmat(NN',1,4),[rsumSC abs(Rt-Rmed) abs(Rt-Rmean) nXN],'.-')
set(gca,'xtick',NN)
xlabel('N')
legend('\Sigma r_{\nu j}','|R_T - R_{median}|','|R_T - R_{mean}|',...
  '||X_N - X_{N-1}||_\infty')



figure(2),clf
plot(repmat(NN',1,3),[Rmax Rmean Rmin],'.-')
hold on
plot(NN([1 end]),repmat(Rt,1,2),'k')
set(gca,'xtick',NN)
xlabel('N')
ylabel('R = log(r_o/r_i) /2\pi')
legend('R_{max}','R_{mean}','R_{min}')
