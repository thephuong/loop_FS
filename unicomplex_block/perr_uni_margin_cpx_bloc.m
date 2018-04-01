function perr = perr_uni_margin_cpx_bloc(m,n,s,rho,rule,SIMUPARAMS)
NTEST = SIMUPARAMS.NTEST;
min_NERR = SIMUPARAMS.min_NERR;
N = m+n;
perr = zeros(1,N-1);
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
for tau=1:N-1
    nerr = 0;
    for ibloc = 1:nbloc
        y = [sbloc; fGenUniVec_cpx_bloc(n,r,ntest_perbloc, -1)] + ...
            1/sqrt(2)*(randn(N,ntest_perbloc)+1i*randn(N,ntest_perbloc));
        %this -tau is tricky to follow convention
        nerr = nerr + sum(fmetric(y,tau,m,n,rho,s) > fmetric(y,0,m,n,rho,s));
        if (nerr > min_NERR)
            break;
        end
    end
    perr(tau) = nerr/ibloc/ntest_perbloc;
end

end
