function func_test_Perr_QAM_subframe(M,N,rho_tot_dB,alpha,NB_SUBFRAME,NTEST)

if mod(N,NB_SUBFRAME) ~= 0
    error('N=%d must divisible to NB_SUBFRAME=%d\n',N,NB_SUBFRAME);
end
if (nargin < 6); NTEST = 1e6; end

fprintf('TEST FS SUBFRAME (QAM=%d): instead of a long SW, use %d short sequences.\n',M,NB_SUBFRAME);

% nopt = 90;mopt = 12;
% N = nopt+mopt;
% rho_tot_dB = -1;
% % rhodB = 0;%0; %3; %rho_d
% alpha = 1; %2; %1;
rho_tot = 10^(rho_tot_dB/10);

epsilon = 1e-3; %1e-6;
stype = 2; %1 uni, 2 ZC, 3 bin
sname = {'uni','ZC','bin'};
% descrip = 'N=102 rho=0db real simu of perr sync, s uni with 3db boost, corr rule. perr margin. epsilon=1e-3';

%CONST
CONST_ML_RULE=1;
CONST_COR_RULE=2;

mm = (3:2:15)*NB_SUBFRAME;
nn = N - mm;
lenmm = length(mm);

ss = cell(lenmm,1);
rhoD_tab = zeros(lenmm,1);

perr = zeros(lenmm,1);         %sync error, optimum rule, sim
ss_sf = cell(lenmm,1);

perr2_iterative1 = zeros(lenmm,1);
% perr2_iterative2 = zeros(lenmm,1);
perr2_iterative4 = zeros(lenmm,1);
perr2_cor_iterative1 = zeros(lenmm,1);
% perr2_cor_iterative2 = zeros(lenmm,1);
perr2_cor_iterative4 = zeros(lenmm,1);

% NB_SUBFRAME=3; %must be impair
DO_VISUALIZE = 0;

k = 6144;

parfor im = 1:lenmm
	m = mm(im);
	n = nn(im);
    msf = m/NB_SUBFRAME;
    nsf = n/NB_SUBFRAME;
    [ss{im},rhoD_tab(im)] = fGenSyncWord(m,n,rho_tot,alpha,stype);
    [ss_sf{im},~] = fGenSyncWord(msf,nsf,rho_tot,alpha,stype);
    
    %real Pe
%     perr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, NTEST);
    perr(im) = perr_QAM_bloc(k,M,m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, NTEST);
    %real Pe iterative ML rule
%     perr2_iterative1(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, NTEST, 1,NB_SUBFRAME);
%     perr2_iterative2(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, NTEST, 2,NBFRAME);
%     perr2_iterative4(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, NTEST, 4,NB_SUBFRAME);
    %real Pe iterative CORR rule
    perr2_cor_iterative1(im) = perr_QAM_iterative_bloc(k,M,msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, NTEST, 1,NB_SUBFRAME);
%     perr2_cor_iterative2(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, NTEST, 2,NBFRAME);
%     perr2_cor_iterative4(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, NTEST, 4,NB_SUBFRAME);
%     fprintf('m=%d perr=%.3e perr2v2=%.3e (trunc=%.3e) : perr2=%.3e\n ===perr2_it_ML  %.3e | %.3e | %.3e\n ===perr2_it_cor %.3e | %.3e | %.3e\n', ...
%         m,perr(im),perr2v2(im),perr2v2_trunc(im),perr2(im), ...
%         perr2_iterative1(im),perr2_iterative2(im),perr2_iterative4(im), ...
%         perr2_cor_iterative1(im),perr2_cor_iterative2(im),perr2_cor_iterative4(im));
    fprintf('SUBFRAME=%d (QAM=%d k=%d) m=%d n=%d --> msf=%d nsf=%d: perr=%.3e\n ===perr2_it_cor %.3e | %.3e\n', ...
        NB_SUBFRAME,M,k,m,n,msf,nsf,perr(im), ...
        perr2_cor_iterative1(im),perr2_cor_iterative4(im));
end

%% Visualize
if (DO_VISUALIZE)
    figure;
    semilogy(mm,perr2v2,'r--','LineWidth',1.5);
    hold on; grid on;
    semilogy(mm,perr,'r:','LineWidth',1.5);
    semilogy(mm,perr2_iterative1,'bo-');
    semilogy(mm,perr2_iterative4,'bd-');
    semilogy(mm,perr2_cor_iterative1,'bo--');
    semilogy(mm,perr2_cor_iterative4,'bd--');
    legend(sprintf('ZC %dm-length',NB_SUBFRAME),'ZC m-length', ...
        'iterative-ML 1 iteration','iterative-ML 4 iterations', ...
        'iterative-Corr 1 iteration','iterative-Corr 4 iterations');
    xlabel('SW length m');
    ylabel('FS error P_e');
end

%% Save
dtt = datestr(datetime);
dtt(dtt==':')='-';
save(sprintf('tpn_SUBFRAME%d_QAM%d_N%d_rho_tot%ddB_alphaPower%d_1e%d_%s.mat',NB_SUBFRAME,M,N,rho_tot_dB,alpha,-log10(epsilon),dtt));