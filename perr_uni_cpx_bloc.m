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
%     switch (rule)
%         case 1
%             metric0 = fmetric_uni_cpx_bloc(y,m,n,rho,s);
%         case 2
%             metric0 = fmetric_cor_cpx_bloc(y,m,s);
%         otherwise
%             metric0 = metric_test(y,m,n,rho,s);
%     end
%     metrictau_bloc = zeros(N-1,ntest_perbloc);
%     switch (rule)
%         case 1
%             for tau = 1:N-1; metrictau_bloc(tau,:) = fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s); end
%         case 2
%             for tau = 1:N-1; metrictau_bloc(tau,:) = fmetric_cor_cpx_bloc(circshift(y,-tau),m,s); end
%         otherwise
%             for tau = 1:N-1; metrictau_bloc(tau,:) = metric_test(circshift(y,-tau),m,n,rho,s); end
%     end
%     metric0_bloc = repmat(metric0,N-1,1);
%     failedtest = sum(metrictau_bloc > metric0_bloc,1) > 0;
%     maybetest = sum(metrictau_bloc == metric0_bloc,1) .* (1-failedtest);
%     maybetest_reduced = maybetest(maybetest > 0);
%     nerr_maybe = 0;
% 	if (~isempty(maybetest_reduced))
%         nerr_maybe = sum(rand(length(maybetest_reduced),1) > 1./maybetest_reduced(:));
%     end
%     nerr = nerr + nerr_maybe + sum(failedtest);

    tau0 = randi(N) - 1;
    y = circshift(y,tau0);
    metrictau_bloc = zeros(N,ntest_perbloc);
    switch (rule)
        case 1
            for tau = 0:N-1; metrictau_bloc(tau+1,:) = fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s); end
        case 2
            for tau = 0:N-1; metrictau_bloc(tau+1,:) = fmetric_cor_cpx_bloc(circshift(y,-tau),m,s); end
        otherwise
            for tau = 0:N-1; metrictau_bloc(tau+1,:) = metric_test(circshift(y,-tau),m,n,rho,s); end
    end
    nerr = nerr + fis_err_bloc(metrictau_bloc,tau0+1);

    if (nerr > 500)
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