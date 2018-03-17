function func_test_Perr_uni_cpx_MULTIFRAME(N,rho_tot_dB,alpha,NBFRAME,NTEST,DO_SAVE)

if nargin < 6; DO_SAVE=1; end

% nopt = 90;mopt = 12;
% N = nopt+mopt;
% rho_tot_dB = -1;
% % rhodB = 0;%0; %3; %rho_d
% alpha = 1; %2; %1;
rho_tot = 10^(rho_tot_dB/10);

epsilon = 1e-3; %1e-6;
stype = 2; %1 uni, 2 ZC, 3 bin
sname = {'uni','ZC','bin'};

%CONST
CONST_ML_RULE=1;
CONST_COR_RULE=2;

SIMUPARAMS.NTEST = NTEST;
SIMUPARAMS.min_NERR = 100;
SIMUPARAMS.CONST_ML_RULE = 1;
SIMUPARAMS.CONST_COR_RULE = 2;

fprintf('TEST FS MULTIFRAMEs: instead of a single frame, use %d FRAMEs.\n', NBFRAME);

mm = 3:2:15; %3:2:17;%floor(N/2);
nn = N - mm;
lenmm = length(mm);

ss = cell(lenmm,1);
rhoD_tab = zeros(lenmm,1);

perr = zeros(lenmm,1);      %sync error, optimum rule, sim
perr2 = zeros(lenmm,1);     %double s
perr2v2 = zeros(lenmm,1);   %one s with double length (+1 because impair ZC)
ssv2 = cell(lenmm,1);
rhoD_tabv2 = zeros(lenmm,1);
% perr2v2_trunc = zeros(lenmm,1); %same as perr2v2 but truncated to have 2m length

perr2_iterative1 = zeros(lenmm,1);
% perr2_iterative2 = zeros(lenmm,1);
perr2_iterative4 = zeros(lenmm,1);
perr2_cor_iterative1 = zeros(lenmm,1);
% perr2_cor_iterative2 = zeros(lenmm,1);
perr2_cor_iterative4 = zeros(lenmm,1);

% NBFRAME=3; %must be impair
DO_VISUALIZE = 0;

% cleanUp function
fname_prefix = sprintf('tpn_MULTIFRAME%d_N%d_rho_tot%ddB_alphaPower%d_1e%d',NBFRAME,N,rho_tot_dB,alpha,-log10(epsilon));
onExitObj = onCleanup(@()onExitFunc(DO_SAVE,fname_prefix));

for im = 1:lenmm
	m = mm(im);
	n = nn(im);
    [ss{im},rhoD_tab(im)] = fGenSyncWord(m,n,rho_tot,alpha,stype);
    if (mod(NBFRAME,2) == 1)
        mv2 = NBFRAME*m;
        nv2 = NBFRAME*n;
    else
        mv2 = NBFRAME*m+1;
        nv2 = NBFRAME*n-1;
    end
    [ssv2{im},rhoD_tabv2(im)] = fGenSyncWord(mv2,nv2,rho_tot,alpha,stype);
    
    %real Pe
    perr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, NTEST);
%     perr2(im) = perr_uni_cpx2_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, NTEST);
    perr2(im) = perr_uni_cpx_bloc(NBFRAME*m,NBFRAME*n,repmat(ss{im},NBFRAME,1),rhoD_tab(im), CONST_ML_RULE, NTEST);
    perr2v2(im) = perr_uni_cpx_bloc(mv2,nv2,ssv2{im},rhoD_tabv2(im), CONST_ML_RULE, NTEST); %%%%%
%     perr2v2_trunc(im) = perr_uni_cpx_bloc(2*m,2*n,ssv2{im}(1:end-1),N/(m*alpha+n)*rho_tot, CONST_ML_RULE, NTEST);
    %real Pe iterative ML rule
    perr2_iterative1(im) = perr_uni_cpx_iterative_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, 1,NBFRAME,SIMUPARAMS);
    perr2_iterative4(im) = perr_uni_cpx_iterative_bloc(m,n,ss{im},rhoD_tab(im), CONST_ML_RULE, 4,NBFRAME,SIMUPARAMS);
    %real Pe iterative CORR rule
    perr2_cor_iterative1(im) = perr_uni_cpx_iterative_bloc(m,n,ss{im},rhoD_tab(im), CONST_COR_RULE, 1,NBFRAME,SIMUPARAMS);
    perr2_cor_iterative4(im) = perr_uni_cpx_iterative_bloc(m,n,ss{im},rhoD_tab(im), CONST_COR_RULE, 4,NBFRAME,SIMUPARAMS);
    fprintf('NBFRAME=%d (randi) m=%d n=%d: perr=%.3e perr2v2=%.3e (mv2=%d nv2=%d) perr2multiple=%.3e\n ===perr2_it_ML  %.3e | %.3e\n ===perr2_it_cor %.3e | %.3e\n', ...
        NBFRAME,m,n,perr(im),perr2v2(im),mv2,nv2,perr2(im), ...
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
    legend(sprintf('ZC %dm-length',NBFRAME),'ZC m-length', ...
        'iterative-ML 1 iteration','iterative-ML 4 iterations', ...
        'iterative-Corr 1 iteration','iterative-Corr 4 iterations');
    xlabel('SW length m');
    ylabel('FS error P_e');
end

%% Save
function onExitFunc(DO_SAVE,fname_prefix)
    if (DO_SAVE)
        dtt = datestr(datetime);
        dtt(dtt==':')='-';
        fname = sprintf('%s_%s.mat',fname_prefix,dtt);
        fprintf('DONE or Interrupted. Saving data to %s\n', fname);
        save(fname);
    end
end

end