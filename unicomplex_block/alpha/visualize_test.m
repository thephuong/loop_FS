clear all;
addpath(genpath(fullfile([fileparts(fileparts(fileparts(pwd))) '/unitfuncs'])));

% % Alpha Uni
% listfile={'tpn_ALPHA_N102_rho_tot0dB_alphaPower1.250000e-01_1e3_14-Mar-2018 01-10-50.mat',...
%     'tpn_ALPHA_N102_rho_tot0dB_alphaPower2.500000e-01_1e3_14-Mar-2018 01-08-51.mat',...
%     'tpn_ALPHA_N102_rho_tot0dB_alphaPower5.000000e-01_1e3_14-Mar-2018 01-13-51.mat',...
% 'tpn_ALPHA_N102_rho_tot0dB_alphaPower1_1e3_13-Mar-2018 20-01-21.mat',...
% 'tpn_N102_rho_tot0dB_alphaPower2_1e3_12-Mar-2018 14-58-06.mat',...
% 'tpn_ALPHA_N102_rho_tot0dB_alphaPower3_1e3_13-Mar-2018 23-35-55.mat',...
% 'tpn_ALPHA_N102_rho_tot0dB_alphaPower4_1e3_13-Mar-2018 23-53-35.mat',...
% 'tpn_ALPHA_N102_rho_tot0dB_alphaPower5_1e3_13-Mar-2018 23-55-50.mat',...
% 'tpn_ALPHA_N102_rho_tot0dB_alphaPower6_1e3_14-Mar-2018 00-00-09.mat'};
%Alpha Normal
listfile={'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower1.250000e-01_1e3_01-Apr-2018 23-27-43.mat', ...
    'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower2.500000e-01_1e3_01-Apr-2018 23-27-50.mat', ...
    'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower5.000000e-01_1e3_01-Apr-2018 23-28-25.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower1_1e3_01-Apr-2018 23-38-27.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower2_1e3_02-Apr-2018 00-48-39.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower3_1e3_02-Apr-2018 02-31-33.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower4_1e3_02-Apr-2018 04-33-30.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower5_1e3_02-Apr-2018 06-44-56.mat', ...
'alpha_Normal/tpn_ALPHA_data_Normal_N102_rho_tot0dB_alphaPower6_1e3_02-Apr-2018 09-05-18.mat'};

nfiles = length(listfile);%size(listfile,1);
%Init
load(listfile{1});
ed_mat = zeros(lenmm,nfiles);
perra_corr_mat = zeros(lenmm,nfiles);
perrcorr_mat = zeros(lenmm,nfiles);
nbits = 51;
debitcorr_tab = zeros(lenmm,nfiles);
debitcorra_tab = zeros(lenmm,nfiles);
debit_ml_tab = zeros(lenmm,nfiles);
debit_ml_a_tab = zeros(lenmm,nfiles);
perr_ml_a_mat = zeros(lenmm,nfiles);
perr_ml_mat = zeros(lenmm,nfiles);
for ifile = 1:nfiles
    load(listfile{ifile});
    %back compatible
    if (~exist('perramargin_ML','var'))
        perr_ml_mat(:,ifile) = perr;
        perr_ml_a_mat(:,ifile) = perr_a_ML;
        perrcorr_mat(:,ifile) = perrcorr;
        perra_corr_mat(:,ifile) = perr_a_cor;
    else
        perr_ml_mat(:,ifile) = perr;
        perr_ml_a_mat(:,ifile) = sum(perramargin_ML,2);
        perrcorr_mat(:,ifile) = perrcorr;
        perra_corr_mat(:,ifile) = sum(perramargin,2);
    end
    
%     k = nbits;
%     nreal = 2*nn(:);
%     rhot = rhoD_tab(:);
%     C = 0.5*log2(1+rhot);
%     V = rhot .* (2+rhot)/2 ./ (1+rhot).^2;
%     ed = qfunc((nreal.*C + 0.5*log2(nreal) - k) ./ sqrt(nreal.*V) / log2(exp(1)));
%     ed_mat=[ed_mat ed];
    ed_mat(:,ifile) = epsilon_D_cpx(nbits,nn,rhoD_tab);
    
    debit_ml_tab(:,ifile) = debit(:);
    debit_ml_a_tab(:,ifile) = debita_ML(:);
    debitcorr_tab(:,ifile) = debitcorr(:);
    debitcorra_tab(:,ifile) = debita(:);
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

if ~exist('data_type_str','var'); data_type_str = 'UniSphere'; end
figure;
semilogy(alpha_tab,pe_ml(im,:),'b-');
hold on; grid on;
semilogy(alpha_tab,pec(im,:),'b--');
semilogy(alpha_tab,pea_ml(im,:),'r*');
semilogy(alpha_tab,pea(im,:),'r+');
xlabel('\alpha = \rho_s/\rho');
ylabel('P_f');
title(sprintf('Frame error P_f m=%d, %s data',mm(im),data_type_str));
legend('Simulation ML rule (f_O)','Estimation ML rule (f_A)','Simulation correlation rule','Estimation correlation rule');