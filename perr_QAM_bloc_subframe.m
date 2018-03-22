%% Generate a long sync word [s1 s2 s3 d1 d2 d3] of length N
% then impose a priori info that frame beginning must occur at somewhere in [s1 d1]
function perr = perr_QAM_bloc_subframe(k,M,m,n,s,rho,rule,NB_SUBFRAME,SIMUPARAMS)

% if nargin < 8; NTEST = 1e6; end
% if nargin < 7; rule = 2; NTEST = 1e6; end

NTEST = SIMUPARAMS.NTEST;
min_NERR = SIMUPARAMS.min_NERR;

rule = SIMUPARAMS.CONST_COR_RULE; %only CORR rule until ML rule is found

if (~ismember(M,[1,2,4,6,8]))
    error('M=%d must be 1,2,4,6,8 !!!');
end

N = m+n;
r = sqrt(n*rho);

switch (rule)
    case SIMUPARAMS.CONST_ML_RULE
        fmetric = @(yy,m,n,rho,s) fmetric_uni_cpx_bloc(yy,m,n,rho,s);
    case SIMUPARAMS.CONST_COR_RULE
        fmetric = @(yy,m,n,rho,s) fmetric_cor_cpx_bloc(yy,m,s);
    otherwise
        error('Rule(%d) not supported.',rule);
end

N_SUBFRAME = N/NB_SUBFRAME;
%a priori information of the position of beginning of subframe
%to be fair in comparison with interative algo
tauval1 = 1;
tauval2 = m/NB_SUBFRAME;
tauval3 = m+1;
tauval4 = m+n/NB_SUBFRAME;
tauvals_tau0 = [tauval1:tauval2 tauval3:tauval4]-1;
EXPECTED_TAUHAT=1;

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);
sbloc = repmat(s,1,ntest_perbloc);
nerr = 0;
for ibloc = 1:nbloc
    y = [sbloc; fGen_lteQAM_Vec_cpx_bloc(k,M,n,r,ntest_perbloc)] + 1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
    tau0 = randi(N) - 1;
    y = circshift(y,tau0);
    %Metric
    tauvals = mod(tauvals_tau0+tau0,N);
    metrictau_bloc = zeros(N_SUBFRAME,ntest_perbloc);
    for itau = 1:length(tauvals)
        tau = tauvals(itau);
        metrictau_bloc(itau,:) = fmetric(circshift(y,-tau),m,n,rho,s);
    end
    %Error calculation
    nerr = nerr + fis_err_bloc(metrictau_bloc,EXPECTED_TAUHAT);
    if (nerr > min_NERR)
        break;
    end
end
perr = nerr/ibloc/ntest_perbloc;

end

function nerr = fis_err_bloc(metrictau,expected_tauhat)
nerr = 0;
for itest = 1:size(metrictau,2)
    nerr = nerr + fis_err(metrictau(:,itest),expected_tauhat);
end
end

function iserr = fis_err(metrictau,expected_tauhat)
maxmetric = max(metrictau);
tauhat_tab = find(metrictau == maxmetric);
tauhat = tauhat_tab(randi(length(tauhat_tab)));
iserr = (tauhat ~= expected_tauhat);
end