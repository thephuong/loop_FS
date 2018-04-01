function perr = perr_uni_margin_cpx_bloc(m,n,s,rho,rule,dataType,SIMUPARAMS)
NTEST = SIMUPARAMS.NTEST;
min_NERR = SIMUPARAMS.min_NERR;
N = m+n;
perr = zeros(1,N-1);
r = sqrt(n*rho);

% Select corresponding Data Generator and Decision metric
[fmetric,fGenData] = anonymous_func_selection(rule,dataType,SIMUPARAMS);

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);

sbloc = repmat(s,1,ntest_perbloc);
for tau=1:N-1
    nerr = 0;
    for ibloc = 1:nbloc
        y = [sbloc; fGenData(n,r,ntest_perbloc)] + ...
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
