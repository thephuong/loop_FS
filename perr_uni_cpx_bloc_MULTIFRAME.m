function perr = perr_uni_cpx_bloc_MULTIFRAME(m,n,s,rho,rule,NBFRAME,SIMUPARAMS)

% if nargin < 6; NTEST = 1e6; end
% if nargin < 5; rule = 2; NTEST = 1e6; end

NTEST = SIMUPARAMS.NTEST;
min_NERR = SIMUPARAMS.min_NERR;

N = m+n;
r = sqrt(n*rho);

switch (rule)
    case SIMUPARAMS.CONST_ML_RULE
        fmetric = @(y,tau,m,n,rho,s) fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s);
    case SIMUPARAMS.CONST_COR_RULE
        fmetric = @(y,tau,m,n,rho,s) fmetric_cor_cpx_bloc(circshift(y,-tau),m,s);
    otherwise
        error('Rule(%d) not supported.',rule);
end

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);
sbloc = repmat(s,1,ntest_perbloc);
nerr = 0;
for ibloc = 1:nbloc
    %Generate data
    x = zeros(NBFRAME*N,ntest_perbloc);
    for iframe = 1:NBFRAME
        x((iframe-1)*N+1:(iframe-1)*N+N,:) = [sbloc;fGenUniVec_cpx_bloc(n,r,ntest_perbloc)];
    end
%     tau0 = randi(nbframe*N) - 1;
    tau0 = 0;
    y = circshift(x + 1/sqrt(2)*(randn(NBFRAME*N,ntest_perbloc)+1i*randn(NBFRAME*N,ntest_perbloc)),tau0);

    %Metric
    tauvals = 0:(N-1); %only if tau0=0
    metrictau_bloc = zeros(N,ntest_perbloc); %get the first FRAME
    for itau = 1:length(tauvals)
        metrictau_bloc(itau,:) = fmetric(y,tauvals(itau),m,n,rho,s);
    end
    nerr = nerr + fis_err_bloc(metrictau_bloc,mod(tau0,N)+1);

    % STOP if enough error
    if (nerr > min_NERR) %500
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
iserr = (~ismember(tauhat,expected_tauhat));
end