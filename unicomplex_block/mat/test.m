load('tpn_N102_rho_tot0dB_alphaPower2_1e3_12-Mar-2018 14-58-06.mat');
% load('tpn_N102_rho_tot-1dB_alphaPower1_1e3_19-Mar-2018 17-48-50.mat');

mmt = mm(:);
nnt = nn(:);

% rhoreal = rhoD_tab; %the power is divised by 2, but noise also then SNR does not change.
% nreal = nnt*2;
% V = rhoreal .* (2+rhoreal)/2 ./ (1+rhoreal).^2;
% C = 0.5*log2(1+rhoreal);

% ep_1e3 = 1e-3;
% R_1e3 = max(0,C - (V./nreal).^0.5*log2(exp(1))*qfuncinv(ep_1e3) + 0.5*log2(nreal)./nreal);
% debitML_1e3 = (1-perr) .* R_1e3 .* nnt;
% debitML_a_1e3 = (1-perraMLc) .* R_1e3 .* nnt;
% debitcorr_1e3 = (1-perrcorr) .* R_1e3 .* nnt;
% debitcorr_a_1e3 = (1-perrac) .* R_1e3 .* nnt;
% 
% nbits = 51;
% ed_51 = qfunc((nnt.*log2(1+rhoD_tab) + 0.5*log2(2*nnt) - nbits) ./ sqrt(nnt.*rhoD_tab.*(rhoD_tab+2)./(rhoD_tab+1).^2));

% Changing name
perrML = perr;
perrML_a = perraMLc;
perrcor = perrcorr;
perrcor_a = perrac;
% debit_epsilon = [Rd,debitML,debitML_a,debitcor,debitcor_a];
% pe_k = [ed,pfML,pfML_a,pfcor,pfcor_a];

epsilon_D = 1.494e-3;
debit_epsilon_D = get_goodput(epsilon_D,nnt,rhoD_tab, perrML,perrML_a,perrcor,perrcor_a);
M = 51;
pe_M = get_pe(M,nnt,rhoD_tab, perrML,perrML_a,perrcor,perrcor_a);
debit_Pf = get_goodput_fixed_Pf(1.806e-2,nnt,rhoD_tab, perrML,perrML_a,perrcor,perrcor_a);

rr = mm/N;
figure;
plot(mm,debit_epsilon_D(:,2),'b-');%,'LineWidth',LWIDTH);
hold on; grid on;
plot(mm,debit_epsilon_D(:,3),'m*');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
plot(mm,debit_epsilon_D(:,4),'b--');%,'LineWidth',LWIDTH);
plot(mm,debit_epsilon_D(:,5),'m+');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
legend('Sim ML','Theory ML','Sim corr','Theory corr');
title(sprintf('epsilon=%.3e N=%d alpha=%d rhoTot=%ddB,s %s,D uni',epsilon_D,N,alpha,rho_tot_dB,sname{stype}));
xlabel('SW length m');
ylabel('G_\epsilon(m,k^*,\epsilon_0)');

figure
subplot(211);
plot(mm,(1-pe_M(:,2))*M,'b-');
hold on; grid on;
plot(mm,(1-pe_M(:,3))*M,'m*');
plot(mm,(1-pe_M(:,4))*M,'b--');
plot(mm,(1-pe_M(:,5))*M,'m+');
axis([xlim 20 55]);
subplot(212);
semilogy(mm,pe_M(:,2),'b-');
hold on; grid on;
semilogy(mm,pe_M(:,3),'m*');
semilogy(mm,pe_M(:,4),'b--');
semilogy(mm,pe_M(:,5),'m+');
% legend('P_f ML','predicted P_f ML','P_f corr','predicted P_f corr');

figure;
semilogy(mm,pe_M(:,2),'b-');
hold on; grid on;
semilogy(mm,pe_M(:,3),'m*');
semilogy(mm,pe_M(:,4),'b--');
semilogy(mm,pe_M(:,5),'m+');
ylabel('P_f');
xlabel('SW length m');
title('Frame error P_f');
legend('Simulation ML rule (f_O)','Estimation ML rule (f_A)','Simulation correlation rule','Estimation correlation rule');

figure;
plot(mm,debit_Pf(:,2),'b-');
hold on; grid on;
plot(mm,debit_Pf(:,3),'m*');
plot(mm,debit_Pf(:,4),'b--');
plot(mm,debit_Pf(:,5),'m+');
legend('Sim ML','Theory ML','Sim corr','Theory corr');
title(sprintf('N=%d alpha=%d rhoTot=%ddB,s %s,D uni',N,alpha,rho_tot_dB,sname{stype}));

figure
semilogy(mm,pe_M(:,1),'bo-');
grid on;
figure
plot(mm,debit_epsilon_D(:,1),'bo-');
grid on;
xlabel('SW length m');
ylabel('k^*(m,\epsilon_0)');
title(sprintf('k^* for epsilon_0=%.3e',epsilon_D))

function debit_Pf = get_goodput_fixed_Pf(Pf,nnt,rhot, ...
    perrML,perrML_a,perrcor,perrcor_a)
    P_FS = 1-Pf;
    epsilon_ML = max(0, 1 - P_FS ./ (1-perrML));
    epsilon_ML_a = max(0, 1 - P_FS ./ (1-perrML_a));
    epsilon_cor = max(0, 1 - P_FS ./ (1-perrcor));
    epsilon_cor_a = max(0, 1 - P_FS ./ (1-perrcor_a));
    
    debit_Pf = zeros(length(nnt),5);
    temp = get_goodput(epsilon_ML,nnt,rhot, perrML,perrML_a,perrcor,perrcor_a);
    debit_Pf(:,2) = temp(:,2);
    temp = get_goodput(epsilon_ML_a,nnt,rhot, perrML,perrML_a,perrcor,perrcor_a);
    debit_Pf(:,3) = temp(:,3);
    temp = get_goodput(epsilon_cor,nnt,rhot, perrML,perrML_a,perrcor,perrcor_a);
    debit_Pf(:,4) = temp(:,4);
    temp = get_goodput(epsilon_cor_a,nnt,rhot, perrML,perrML_a,perrcor,perrcor_a);
    debit_Pf(:,5) = temp(:,5);
end

function debit_epsilon = get_goodput(epsilon,nnt,rhot, ...
    perrML,perrML_a,perrcor,perrcor_a)
    nreal = 2*nnt(:);
    V = rhot .* (2+rhot)/2 ./ (1+rhot).^2;
    C = 0.5*log2(1+rhot);
    Md = max(0, nreal .* C - ...
        sqrt(V.*nreal)*log2(exp(1)).*qfuncinv(epsilon) + ...
        0.5*log2(nreal));
    debitML = (1-perrML) .* Md;
    debitML_a = (1-perrML_a) .* Md;
    debitcor = (1-perrcor) .* Md;
    debitcor_a = (1-perrcor_a) .* Md;
    debit_epsilon = [Md,debitML,debitML_a,debitcor,debitcor_a];
end

function pe_k = get_pe(k,nnt,rhot, perrML,perrML_a,perrcor,perrcor_a)
    nreal = 2*nnt(:);
    C = 0.5*log2(1+rhot);
    V = rhot .* (2+rhot)/2 ./ (1+rhot).^2;
    ed = qfunc((nreal.*C + 0.5*log2(nreal) - k) ./ sqrt(nreal.*V) / log2(exp(1)));
    pfML = 1-(1-perrML).*(1-ed);
    pfML_a = 1-(1-perrML_a).*(1-ed);
    pfcor = 1-(1-perrcor).*(1-ed);
    pfcor_a = 1-(1-perrcor_a).*(1-ed);
    pe_k = [ed,pfML,pfML_a,pfcor,pfcor_a];
end