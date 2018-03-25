clear all;
listfile={'tpn_ALPHA_N102_rho_tot0dB_alphaPower1.250000e-01_1e3_14-Mar-2018 01-10-50.mat',...
    'tpn_ALPHA_N102_rho_tot0dB_alphaPower2.500000e-01_1e3_14-Mar-2018 01-08-51.mat',...
    'tpn_ALPHA_N102_rho_tot0dB_alphaPower5.000000e-01_1e3_14-Mar-2018 01-13-51.mat',...
'tpn_ALPHA_N102_rho_tot0dB_alphaPower1_1e3_13-Mar-2018 20-01-21.mat',...
'tpn_N102_rho_tot0dB_alphaPower2_1e3_12-Mar-2018 14-58-06.mat',...
'tpn_ALPHA_N102_rho_tot0dB_alphaPower3_1e3_13-Mar-2018 23-35-55.mat',...
'tpn_ALPHA_N102_rho_tot0dB_alphaPower4_1e3_13-Mar-2018 23-53-35.mat',...
'tpn_ALPHA_N102_rho_tot0dB_alphaPower5_1e3_13-Mar-2018 23-55-50.mat',...
'tpn_ALPHA_N102_rho_tot0dB_alphaPower6_1e3_14-Mar-2018 00-00-09.mat'};
nfiles = length(listfile);%size(listfile,1);
ed_mat = [];
perra_corr_mat = [];
perrcorr_mat = [];
nbits = 51;
debitcorr_tab = [];
debitcorra_tab = [];
debit_ml_tab = [];
debit_ml_a_tab = [];
perr_ml_a_mat = [];
perr_ml_mat = [];
for ifile = 1:nfiles
    load(listfile{ifile});
    perr_ml_mat = [perr_ml_mat perr];
    perr_ml_a_mat = [perr_ml_a_mat sum(perramargin_ML,2)];
    
    perrcorr_mat=[perrcorr_mat perrcorr];
    perra_corr_mat=[perra_corr_mat sum(perramargin,2)];
%     epsilonD = qfunc((nn(:).*log2(1+rhoD_tab(:)) + 0.5*log2(2*nn(:)) - nbits) ./ sqrt(nn(:).*rhoD_tab(:).*(rhoD_tab(:)+2)./(rhoD_tab(:)+1).^2));
    ed_mat=[ed_mat ed];
    
    debit_ml_tab = [debit_ml_tab debit(:)];
    debit_ml_a_tab = [debit_ml_a_tab debita_ML(:)];
    debitcorr_tab = [debitcorr_tab debitcorr(:)];
    debitcorra_tab = [debitcorra_tab debita(:)];
end
pea = 1-(1-min(perra_corr_mat,1)).*(1-ed_mat);
pec = 1-(1-min(perrcorr_mat,1)).*(1-ed_mat);

pea_ml = 1-(1-min(perr_ml_a_mat,1)).*(1-ed_mat);
pe_ml = 1-(1-min(perr_ml_mat,1)).*(1-ed_mat);

% figure;
% semilogy(mm,pec);
% figure;
% semilogy(mm,pea);
% figure;
% semilogy([0.5 1:8],min(pea));

% figure;
% plot(mm,debitcorr_tab);
% figure;
% plot(mm,debitcorra_tab);
figure;
plot(debitcorr_tab.');
figure;
plot(debitcorra_tab.');

alpha_tab = [0.125 0.25 0.5 1:6];
im = 7;
figure;
subplot(211);
plot(alpha_tab,debit_ml_tab(im,:),'b-');
hold on; grid on;
plot(alpha_tab,debitcorr_tab(im,:),'b--');
plot(alpha_tab,debit_ml_a_tab(im,:),'r*');
plot(alpha_tab,debitcorra_tab(im,:),'r+');
xlabel('\alpha = \rho_s/\rho');
title('(a) Goodput G');
ylabel('bits/frame');

subplot(212);
semilogy(alpha_tab,pe_ml(im,:),'b-');
hold on; grid on;
semilogy(alpha_tab,pec(im,:),'b--');
semilogy(alpha_tab,pea_ml(im,:),'r*');
semilogy(alpha_tab,pea(im,:),'r+');
xlabel('\alpha = \rho_s/\rho');
title('(b) Overall frame error P_f');