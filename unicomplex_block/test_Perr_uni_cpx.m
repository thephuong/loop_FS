clear all
addpath(genpath('../unitfuncs/'));

nopt = 90;mopt = 12;
N = nopt+mopt;
rho_tot_dB = -1;
IS_PREDICT_ONLY = 3;
% rhodB = 0;%0; %3; %rho_d
alpha = 1; %2; %1;
rho_tot = 10^(rho_tot_dB/10);

epsilon = 1e-3; %1e-6;
stype = 2; %1 uni, 2 ZC, 3 bin
sname = {'uni','ZC','bin'};
% descrip = 'N=102 rho=0db real simu of perr sync, s uni with 3db boost, corr rule. perr margin. epsilon=1e-3';

CONST_ML_RULE=1;
CONST_COR_RULE=2;
CONST_dataType_UNISPHERE=1;
CONST_dataType_GAUSSIAN=2;
dataTypeName = {'UniSphere','Normal'};

NTEST = 1e5;
mm = 3:2:28;%3:2:32;%floor(N/2);
nn = N - mm;
lenmm = length(mm);
DATA_TYPE=CONST_dataType_UNISPHERE;

SIMUPARAMS = struct('NTEST',NTEST, ...
    'CONST_ML_RULE',CONST_ML_RULE,'CONST_COR_RULE',CONST_COR_RULE, ...
    'CONST_dataType_UNISPHERE',CONST_dataType_UNISPHERE,'CONST_dataType_GAUSSIAN',CONST_dataType_GAUSSIAN, ...
    'min_NERR', 200);

if (stype == 2)
    pm = find(mod(mm,2)==1);
    mm = mm(pm);
    nn = N-mm;
end
ss = cell(lenmm,1);
rhoD_tab = zeros(lenmm,1);

% perr_MLtest = zeros(length(mm),1);         %sync error, optimum rule, sim
perr = zeros(lenmm,1);         %sync error, optimum rule, sim
perrcorr = zeros(lenmm,1);     %sync error, corr rule, sim

perr_margin_ML = zeros(lenmm,N-1);   %Pe at tau, sim, ML
perr_margin_cor = zeros(lenmm,N-1);      %Pe at tau, sim, corr rule

perr_a_cor = zeros(lenmm,1);            %Pe,union, theory, corr rule
perr_a_ML = zeros(lenmm,1);             %Pe,union, theory, ML rule
perr_a_margin_cor = zeros(lenmm,N-1);%Pe at tau, theory
perr_a_margin_ML = zeros(lenmm,N-1);%Pe at tau, theory, ML

parfor im = 1:lenmm
	m = mm(im);
	n = nn(im);
    fprintf('m=%d\n',m);
    [ss{im},rhoD_tab(im)] = fGenSyncWord(m,n,rho_tot,alpha,stype);
    
    %real frame sync (FS) error
    if (IS_PREDICT_ONLY >= 2)
%         perr(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 1, NTEST);
%         perrcorr(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 2, NTEST);
%     %     perr_MLtest(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 3, NTEST);
        perr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_ML_RULE,DATA_TYPE, SIMUPARAMS);
        perrcorr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_COR_RULE,DATA_TYPE, SIMUPARAMS);
    end
    
    %real Pe(tau)
    if (IS_PREDICT_ONLY >= 3)
%         perrmargin_ML(im,:) = perr_uni_margin_ML_cpx(m,n,ss{im},rhoD_tab(im),NTEST);
%         perrmargin(im,:) = perr_uni_margin_corr_cpx(m,n,ss{im},rhoD_tab(im),NTEST);
        perr_margin_ML(im,:) = perr_uni_margin_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_ML_RULE,DATA_TYPE, SIMUPARAMS);
        perr_margin_cor(im,:) = perr_uni_margin_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_COR_RULE,DATA_TYPE, SIMUPARAMS);
    end

    %theoretical Pe(tau)
    if (IS_PREDICT_ONLY >= 1)
        if (DATA_TYPE == CONST_dataType_UNISPHERE)
            perr_a_margin_cor(im,:) = fpredict_perr_uni_margin_corr_cpx(ss{im},m,n,rhoD_tab(im));
            perr_a_cor(im) = sum(perr_a_margin_cor(im,:),2);
            perr_a_margin_ML(im,:) = fpredict_perr_uni_margin_ML_cpx(ss{im},m,n,rhoD_tab(im));
            perr_a_ML(im) = sum(perr_a_margin_ML(im,:),2);
        else
            perr_a_margin_cor(im,:) = fpredict_perr_Normal_margin_corr_cpx(ss{im},m,n,rhoD_tab(im));
            perr_a_cor(im) = sum(perr_a_margin_cor(im,:),2);
            perr_a_margin_ML(im,:) = fpredict_perr_Normal_margin_ML_cpx(ss{im},m,n,rhoD_tab(im));
            perr_a_ML(im) = sum(perr_a_margin_ML(im,:),2);
        end
    end

	fprintf('m=%d DONE ..\n',m);
end

nnt = nn(:);
R = k_D_cpx(epsilon,nn,rhoD_tab) ./ nnt;
nbits = 51;
ed = epsilon_D_cpx(nbits,nn,rhoD_tab);

% %% Data
% rhoreal = rho/2;
% nnreal = nn*2;
% 
% V = rhoreal*(2+rhoreal)/2/(1+rhoreal)^2;
% C = 0.5*log2(1+rhoreal);
% R = C - (V./nnreal).^0.5*log2(exp(1))*qfuncinv(epsilon) + 0.5*log2(nnreal)./nnreal;

%% Visualize
Rc = max(0,R);
perrac = min(1,perr_a_cor);
perraMLc = min(1,perr_a_ML);
debit = (1-perr(:)) .* Rc .* nn(:);
debitcorr = (1-perrcorr(:)) .* Rc .* nnt;
debita = (1-perrac(:)) .* Rc .* nnt;
debita_ML = (1-perraMLc(:)) .* Rc .* nnt;

MKERSIZE = 5; MKERSIZE_BIG = 5;
LWIDTH = 1; LWIDTH_BOLD = 1;
figure;
semilogy(mm,perr,'b-');%,'LineWidth',LWIDTH);
hold on; grid on;
semilogy(mm,sum(perr_margin_ML,2),'ro');%,'LineWidth',LWIDTH,'MarkerSize',MKERSIZE_BIG);
semilogy(mm,perr_a_ML,'m*');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
semilogy(mm,perrcorr,'b--','LineWidth', LWIDTH);
semilogy(mm,sum(perr_margin_cor,2),'rs');%,'LineWidth',LWIDTH,'MarkerSize',MKERSIZE_BIG);
semilogy(mm,perr_a_cor,'m+');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
legend('Pe ML rule','Pe,union ML (sim)','Pe,union ML (predicted)',...
    'Pe cor rule','Pe,union cor (sim)','Pe,union cor (predicted)');
title(sprintf('N=%d alpha=%d rhoTot=%ddB,s %s,D uni',N,alpha,rho_tot_dB,sname{stype}));

figure;
plot(mm,debit,'b-');%,'LineWidth',LWIDTH);
hold on; grid on;
plot(mm,debita_ML,'m*');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
plot(mm,debitcorr,'b--');%,'LineWidth',LWIDTH);
plot(mm,debita,'m+');%,'LineWidth',LWIDTH_BOLD,'MarkerSize',MKERSIZE);
legend('Sim ML','Theory ML','Sim corr','Theory corr');
title(sprintf('N=%d alpha=%d rhoTot=%ddB,s %s,D uni',N,alpha,rho_tot_dB,sname{stype}));

figure
semilogy(mm,1-(1-perr).*(1-ed));
hold on; grid on;
semilogy(mm,1-(1-min(perr_a_ML,1)).*(1-ed),'m*');
semilogy(mm,1-(1-perrcorr).*(1-ed),'b--');
semilogy(mm,1-(1-min(perr_a_cor,1)).*(1-ed),'m+');
legend('P_f ML','predicted P_f ML','P_f corr','predicted P_f corr');

%Corr rule
figure; im=3;
semilogy(perr_margin_cor(im,:),'bo-');
hold on; grid on;
semilogy(perr_a_margin_cor(im,:),'r*--');
title(['Corr rule. P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');
%Corr rule
figure; im=10;
semilogy(perr_margin_cor(im,:),'bo-');
hold on; grid on;
semilogy(perr_a_margin_cor(im,:),'r*--');
title(['Corr rule. P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');

%ML rule
figure; im=3;
semilogy(perr_margin_ML(im,:),'bo-');
hold on; grid on;
semilogy(perr_a_margin_ML(im,:),'r*--');
title(['ML rule. P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');
%ML rule
figure; im=10;
semilogy(perr_margin_ML(im,:),'bo-');
hold on; grid on;
semilogy(perr_a_margin_ML(im,:),'r*--');
title(['ML rule. P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');

%% Save
if (IS_PREDICT_ONLY >= 3)
    dtt = datestr(datetime);
    dtt(dtt==':')='-';
    save(sprintf('tpn_data_%s_N%d_rho_tot%ddB_alphaPower%d_1e%d_%s.mat',dataTypeName{DATA_TYPE}, ...
        N,rho_tot_dB,alpha,-log10(epsilon),dtt));
end