function perr = perr_uni_margin_corr_cpx(m,n,s,rho,NTEST,modeSIM)
if nargin < 4; NTEST=1e6; end
if nargin < 6; modeSIM=1; end %while mode sequential
N = m+n;
perr = zeros(1,N-1);
r = sqrt(n*rho);
if (modeSIM == 1)
    for tau=1:N-1
        nerr = 0;
        for i = 1:NTEST
            y = [s; fGenUniVec_cpx(n,r)] + 1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
            yy = circshift(y,-tau); %this -tau is tricky to follow convention
            nerr = nerr + (real(yy(1:m)' * s) > real(y(1:m)' * s));
            if (nerr > 300); break; end
        end
        perr(tau) = nerr/i;
    end
else
    parfor tau=1:N-1
        nerr = zeros(NTEST,1);
        for i = 1:NTEST
            y = [s; fGenUniVec_cpx(n,r)] + 1/sqrt(2)*(randn(N,1)+1i*randn(N,1));
            yy = circshift(y,-tau); %this -tau is tricky to follow convention
            nerr(i) = (real(yy(1:m)'*s) > real(y(1:m)'*s));
        end
        perr(tau) = mean(nerr);
    end
end
