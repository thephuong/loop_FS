function perr = perr_uni_cpx_bloc(m,n,s,rho,rule,dataType,SIMUPARAMS)

% if nargin < 6; NTEST = 1e6; end
% if nargin < 5; rule = 2; NTEST = 1e6; end
NTEST = SIMUPARAMS.NTEST;
min_NERR = SIMUPARAMS.min_NERR;

N = m+n;
r = sqrt(n*rho);

switch (dataType)
    case SIMUPARAMS.CONST_dataType_UNISPHERE
        fGenData = @(n,r,ntest_perbloc) fGenUniVec_cpx_bloc(n,r,ntest_perbloc,-1);
        switch (rule)
            case SIMUPARAMS.CONST_ML_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s);
            case SIMUPARAMS.CONST_COR_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_cor_cpx_bloc(circshift(y,-tau),m,s);
            otherwise
                error('Rule(%d) not supported.',rule);
        end
    case SIMUPARAMS.CONST_dataType_GAUSSIAN
        fGenData = @(n,r,ntest_perbloc) r/sqrt(2*n)*(randn(n,ntest_perbloc) + 1i*randn(n,ntest_perbloc));
        switch (rule)
            case SIMUPARAMS.CONST_ML_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_Normal_cpx_bloc(circshift(y,-tau),m,rho,s);
            case SIMUPARAMS.CONST_COR_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_cor_cpx_bloc(circshift(y,-tau),m,s);
            otherwise
                error('Rule(%d) not supported.',rule);
        end
    otherwise
        error('dataType(%d) not supported.',dataType);
end

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);
sbloc = repmat(s,1,ntest_perbloc);
nerr = 0;
for ibloc = 1:nbloc
    y = [sbloc; fGenData(n,r,ntest_perbloc)] + ...
        1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
    metric0_bloc = repmat(fmetric(y,0,m,n,rho,s),N-1,1);
    metrictau_bloc = zeros(N-1,ntest_perbloc);
    for tau = 1:N-1
        metrictau_bloc(tau,:) = fmetric(y,tau,m,n,rho,s);
    end
    failedtest = sum(metrictau_bloc > metric0_bloc,1) > 0;
    maybetest = (sum(metrictau_bloc == metric0_bloc) > 0) .* (1-failedtest);
    nerr = nerr + sum(rand(sum(maybetest),1) > 0.5) + sum(failedtest);
    if (nerr > min_NERR)
        break;
    end
end
perr = nerr/ibloc/ntest_perbloc;

end

function metric = metric_test(y,m,n,rho,s)
ys = y(1:m);
A = sqrt(n*rho)/norm(y) - (2*n+1)/4/norm(y)^2;
metric = -norm(ys - s/A)^2;
end

% function metric = metric_uni(y,m,n,rho,s)
% ys = y(1:m);
% yd = y(m+1:end);
% metric = 2*real(ys'*s) + log(besseli(n-1,2*sqrt(n*rho)*norm(yd))) - n*log(norm(yd));
% end

% function metric = metric_cor(y,m,s)
% ys = y(1:m);
% metric = real(ys' * s);
% end