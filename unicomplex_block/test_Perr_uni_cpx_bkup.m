clear all
nopt = 90;mopt = 12;
N = nopt+mopt;
rho_tot_dB = 0;
% rhodB = 0;%0; %3; %rho_d
alpha = 1; %2; %1;
rho_tot = 10^(rho_tot_dB/10);

epsilon = 1e-3; %1e-6;
stype = 2; %1 uni, 2 ZC, 3 bin
sname = {'uni','ZC','bin'};
% descrip = 'N=102 rho=0db real simu of perr sync, s uni with 3db boost, corr rule. perr margin. epsilon=1e-3';

CONST_ML_RULE=1;
CONST_COR_RULE=2;

NTEST = 1e6;
mm = 3:2:28;%floor(N/2);
nn = N - mm;
lenmm = length(mm);

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

perrmargin_ML = zeros(lenmm,N-1);   %Pe at tau, sim, ML
perrmargin = zeros(lenmm,N-1);      %Pe at tau, sim, corr rule

perra = zeros(lenmm,1);        %Pe,union, theory
perramargin = zeros(lenmm,N-1);%Pe at tau, theory
perramargin_ML = zeros(lenmm,N-1);%Pe at tau, theory, ML
R = zeros(1,lenmm);

parfor im = 1:lenmm
	m = mm(im);
	n = nn(im);
    fprintf('m=%d\n',m);
    [ss{im},rhoD_tab(im)] = fGenSyncWord(m,n,rho_tot,alpha,stype);
    
    %real Pe
    perr(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 1, NTEST);
    perrcorr(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 2, NTEST);
% %     perr_MLtest(im) = perr_uni_cpx(m,n,ss{im},rhoD_tab(im), 3, NTEST);
%     perr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_ML_RULE, NTEST);
%     perrcorr(im) = perr_uni_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_COR_RULE, NTEST);
    
    %real Pe(tau)
%     perrmargin_ML(im,:) = perr_uni_margin_ML_cpx(m,n,ss{im},rhoD_tab(im),NTEST);
%     perrmargin(im,:) = perr_uni_margin_corr_cpx(m,n,ss{im},rhoD_tab(im),NTEST);
    perrmargin_ML(im,:) = perr_uni_margin_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_ML_RULE,NTEST);
    perrmargin(im,:) = perr_uni_margin_cpx_bloc(m,n,ss{im},rhoD_tab(im),CONST_COR_RULE,NTEST);
    
    %theoretical Pe(tau)
    perramargin(im,:) = fpredict_perr_uni_margin_corr_cpx(ss{im},m,n,rhoD_tab(im));
    perra(im) = sum(perramargin(im,:),2);
    perramargin_ML(im,:) = fpredict_perr_uni_margin_ML_cpx(ss{im},m,n,rhoD_tab(im));
    
%     rhoreal = rho/2;
    rhoreal = rhoD_tab(im); %the power is divised by 2, but noise also then SNR does not change.
    nreal = n*2;
    V = rhoreal*(2+rhoreal)/2/(1+rhoreal)^2;
    C = 0.5*log2(1+rhoreal);
    R(im) = C - (V./nreal).^0.5*log2(exp(1))*qfuncinv(epsilon) + 0.5*log2(nreal)./nreal;
end

%nbits = 16;
%ed = qfunc((nn.*log2(1+rhoD_tab) + 0.5*log2(2*nn) - nbits) ./ sqrt(nn.*rhoD_tab.*(rhoD_tab+2)./(rhoD_tab+1).^2));

% %% Data
% rhoreal = rho/2;
% nnreal = nn*2;
% 
% V = rhoreal*(2+rhoreal)/2/(1+rhoreal)^2;
% C = 0.5*log2(1+rhoreal);
% R = C - (V./nnreal).^0.5*log2(exp(1))*qfuncinv(epsilon) + 0.5*log2(nnreal)./nnreal;

%% Visualize
Rc = max(0,R);
perrac = min(1,perra);
debit = (1-perr).' .* Rc .* nn;
debitcorr = (1-perrcorr).' .* Rc .* nn;
debita = (1-perrac).' .* Rc .* nn;

figure;
semilogy(mm,perrcorr,'b-');
hold on; grid on;
% semilogy(mm,perr,'b--');
semilogy(mm,sum(perrmargin,2),'ro');
semilogy(mm,perra,'m*');
legend('Sync err, corr rule','Pe,union sim','Pe,union theory');
title(sprintf('N=%d alpha=%d rhoTot=%ddB,s %s,D uni',N,alpha,rho_tot_dB,sname{stype}));

figure;
plot(mm,debitcorr,'b-');
hold on; grid on;
plot(mm,debita,'m*');
legend('Sim','Theory');
title(sprintf('N=%d alpha=%d rhoTot=%ddB,s %s,D uni',N,alpha,rho_tot_dB,sname{stype}));

figure; im=3;
semilogy(perrmargin(im,:),'bo-');
hold on; grid on;
semilogy(perramargin(im,:),'r*--');
title(['P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');

figure; im=10;
semilogy(perrmargin(im,:),'bo-');
hold on; grid on;
semilogy(perramargin(im,:),'r*--');
title(['P_e(\tau).' sprintf('m=%d n=%d alpha=%d rhoTotdB=%d',mm(im),nn(im),alpha,rho_tot_dB)]);
legend('Sim','Theory');

%% Save
dtt = datestr(datetime);
dtt(dtt==':')='-';
save(sprintf('tpn_N%d_rho_tot%ddB_alphaPower%d_1e%d_%s.mat',N,rho_tot_dB,alpha,-log10(epsilon),dtt));