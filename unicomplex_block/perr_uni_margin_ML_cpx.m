function perr = perr_uni_margin_ML_cpx(m,n,s,rho,NTEST)
if nargin < 4; NTEST=1e6; end
N = m+n;
perr = zeros(1,N-1);
r = sqrt(n*rho);

ntest_perbloc = 100;
nbloc = ceil(NTEST/ntest_perbloc);

sbloc = repmat(s,1,ntest_perbloc);
for tau=1:N-1
    nerr = 0;
    for i = 1:nbloc
        y = [sbloc; fGenUniVec_cpx(n,r,ntest_perbloc)] + 1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
        %this -tau is tricky to follow convention
        nerr = nerr + sum(metric_uni(circshift(y,-tau),m,n,rho,s) > metric_uni(y,m,n,rho,s));
        if (nerr > 50); break; end
    end
    perr(tau) = nerr/i/ntest_perbloc;
end

end

function metric = metric_uni(y,m,n,rho,s)
% ys = y(1:m,:);
normyd = sqrt(sum(y(m+1:end,:).^2,1));
metric = 2*real(s'*y(1:m,:)) + log(besseli(n-1,2*sqrt(n*rho)*normyd)) - n*log(normyd);
end
