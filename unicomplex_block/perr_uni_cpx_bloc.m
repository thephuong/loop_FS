function perr = perr_uni_cpx_bloc(m,n,s,rho,rule,NTEST)

if nargin < 6; NTEST = 1e6; end
if nargin < 5; rule = 2; NTEST = 1e6; end
N = m+n;
r = sqrt(n*rho);

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);
sbloc = repmat(s,1,ntest_perbloc);
nerr = 0;
for ibloc = 1:nbloc
    y = [sbloc; fGenUniVec_cpx_bloc(n,r,ntest_perbloc)] + 1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
    switch (rule)
        case 1
            metric0 = fmetric_uni_cpx_bloc(y,m,n,rho,s);
        case 2
            metric0 = fmetric_cor_cpx_bloc(y,m,s);
        otherwise
            metric0 = metric_test(y,m,n,rho,s);
    end
    metric0_bloc = repmat(metric0,N-1,1);
    metrictau_bloc = zeros(N-1,ntest_perbloc);
    switch (rule)
        case 1
            for tau = 1:N-1; metrictau_bloc(tau,:) = fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s); end
        case 2
            for tau = 1:N-1; metrictau_bloc(tau,:) = fmetric_cor_cpx_bloc(circshift(y,-tau),m,s); end
        otherwise
            for tau = 1:N-1; metrictau_bloc(tau,:) = metric_test(circshift(y,-tau),m,n,rho,s); end
    end
    failedtest = sum(metrictau_bloc > metric0_bloc,1) > 0;
    maybetest = (sum(metrictau_bloc == metric0_bloc) > 0) .* (1-failedtest);
    nerr = nerr + sum(rand(sum(maybetest),1) > 0.5) + sum(failedtest);
    if (nerr > 200)
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