function perr = perr_uni_margin_cpx_bloc(m,n,s,rho,rule,NTEST)
if nargin < 4; NTEST=1e6; end
N = m+n;
perr = zeros(1,N-1);
r = sqrt(n*rho);

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);

sbloc = repmat(s,1,ntest_perbloc);
for tau=1:N-1
    nerr = 0;
    for ibloc = 1:nbloc
        y = [sbloc; fGenUniVec_cpx_bloc(n,r,ntest_perbloc)] + 1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
        %this -tau is tricky to follow convention
        switch (rule)
            case 1
                nerr = nerr + sum(fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s) > fmetric_uni_cpx_bloc(y,m,n,rho,s));
            otherwise
                nerr = nerr + sum(fmetric_cor_cpx_bloc(circshift(y,-tau),m,s) > fmetric_cor_cpx_bloc(y,m,s));
        end
        if (nerr > 200)
            break;
        end
    end
    perr(tau) = nerr/ibloc/ntest_perbloc;
end

end
