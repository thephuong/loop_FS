function func_test_Perr_uni_cpx_subframe(N,rho_tot_dB,alpha,NB_SUBFRAME,NTEST)

if mod(N,NB_SUBFRAME) ~= 0
    error('N=%d must divisible to NB_SUBFRAME=%d\n',N,NB_SUBFRAME);
end
if (nargin < 5); NTEST = 1e6; end

%CONST
CONST_ML_RULE=1;
CONST_COR_RULE=2;

SIMUPARAMS.NTEST = NTEST;
SIMUPARAMS.min_NERR = 100;
SIMUPARAMS.CONST_ML_RULE = 1;
SIMUPARAMS.CONST_COR_RULE = 2;

fprintf('TEST FS SUBFRAME: instead of a long SW, use %d short sequences.\n', NB_SUBFRAME);

% nopt = 90;mopt = 12;
% N = nopt+mopt;
% rho_tot_dB = -1;
% % rhodB = 0;%0; %3; %rho_d
% alpha = 1; %2; %1;
rho_tot = 10^(rho_tot_dB/10);

epsilon = 1e-3; %1e-6;
stype = 2; %1 uni, 2 ZC, 3 bin
sname = {'uni','ZC','bin'};

% cleanUp function
fname_prefix = sprintf('tpn_SUBFRAME%d_N%d_rho_tot%ddB_alphaPower%d_1e%d',NB_SUBFRAME,N,rho_tot_dB,alpha,-log10(epsilon));
onExitObj = onCleanup(@()onExitFunc(fname_prefix));

mm = (3:2:15)*NB_SUBFRAME;
nn = N - mm;
lenmm = length(mm);

ss = cell(lenmm,1);
rhoD_tab = zeros(lenmm,1);

perr = zeros(lenmm,1);         %sync error, optimum rule, sim
perr_cor = zeros(lenmm,1);         %sync error, optimum rule, sim
ss_sf = cell(lenmm,1);

perr2_iterative1 = zeros(lenmm,1);
% perr2_iterative2 = zeros(lenmm,1);
perr2_iterative4 = zeros(lenmm,1);
perr2_cor_iterative1 = zeros(lenmm,1);
% perr2_cor_iterative2 = zeros(lenmm,1);
perr2_cor_iterative4 = zeros(lenmm,1);

% NB_SUBFRAME=3; %must be impair
DO_VISUALIZE = 0;

for im = 1:lenmm
	m = mm(im);
	n = nn(im);
    msf = m/NB_SUBFRAME;
    nsf = n/NB_SUBFRAME;
    [ss{im},rhoD_tab(im)] = fGenSyncWord(m,n,rho_tot,alpha,stype);
    [ss_sf{im},~] = fGenSyncWord(msf,nsf,rho_tot,alpha,stype);
    
    %real Pe
%     perr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, NTEST);
    perr(im) = perr_uni_cpx_bloc_subframe(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE,NB_SUBFRAME,SIMUPARAMS);
    perr_cor(im) = perr_uni_cpx_bloc_subframe(m,n,ss{im},rhoD_tab(im), CONST_COR_RULE,NB_SUBFRAME,SIMUPARAMS);
    %real Pe iterative ML rule
    perr2_iterative1(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, 1,NB_SUBFRAME,SIMUPARAMS);
%     perr2_iterative2(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, 2,NB_SUBFRAME,SIMUPARAMS);
    perr2_iterative4(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_ML_RULE, 4,NB_SUBFRAME,SIMUPARAMS);
    %real Pe iterative CORR rule
    perr2_cor_iterative1(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, 1,NB_SUBFRAME,SIMUPARAMS);
%     perr2_cor_iterative2(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, 2,NB_SUBFRAME,SIMUPARAMS);
    perr2_cor_iterative4(im) = perr_uni_cpx_iterative_bloc(msf,nsf,ss_sf{im},rhoD_tab(im), CONST_COR_RULE, 4,NB_SUBFRAME,SIMUPARAMS);
    fprintf('SUBFRAME=%d (randi) m=%d n=%d --> msf=%d nsf=%d: perr=%.3e perr_cor=%.3e\n ===perr2_it_ML  %.3e | %.3e\n ===perr2_it_cor %.3e | %.3e\n', ...
        NB_SUBFRAME,m,n,msf,nsf,perr(im),perr_cor(im), ...
        perr2_iterative1(im),perr2_iterative4(im), ...
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
function onExitFunc(fname_prefix)
    dtt = datestr(datetime);
    dtt(dtt==':')='-';
    fname = sprintf('%s_%s.mat',fname_prefix,dtt);
    fprintf('DONE or Interrupted. Saving data to %s\n', fname);
    save(fname);
end

end