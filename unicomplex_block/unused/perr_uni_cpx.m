function perr = perr_uni_cpx(m,n,s,rho,rule,NTEST)

if nargin < 6; NTEST = 1e6; end
if nargin < 5; rule = 2; NTEST = 1e6; end
N = m+n;
r = sqrt(n*rho);

nerr = zeros(NTEST,1);
% nerrs = 0;
for i = 1:NTEST
    y = [s; fGenUniVec_cpx(n,r)] + 1/sqrt(2)*(randn(N,1)+1i*randn(N,1));

    switch (rule)
        case 1
            metric0 = metric_uni(y,m,n,rho,s);
        case 2
            metric0 = metric_cor(y,m,s);
        otherwise
            metric0 = metric_test(y,m,n,rho,s);
    end
    posm = 0;
    for tau = 1:N-1
        switch (rule)
            case 1
                metrictau = metric_uni(circshift(y,-tau),m,n,rho,s);
            case 2
                metrictau = metric_cor(circshift(y,-tau),m,s);
            otherwise
                metrictau = metric_test(circshift(y,-tau),m,n,rho,s);
        end
        if (metrictau > metric0)
            nerr(i) = 1;
            break;
        elseif (metrictau == metric0)
            posm = tau;
        end
    end
    
    if (posm > 0); nerr(i) = (rand() > 0.5); end
    if (sum(nerr) > 300)
        break;
    end
end

% perr = mean(nerr);
perr = sum(nerr)/i;
end

function metric = metric_test(y,m,n,rho,s)
ys = y(1:m);
A = sqrt(n*rho)/norm(y) - (2*n+1)/4/norm(y)^2;
metric = -norm(ys - s/A)^2;
end

function metric = metric_uni(y,m,n,rho,s)
ys = y(1:m);
yd = y(m+1:end);
metric = 2*real(ys'*s) + log(besseli(n-1,2*sqrt(n*rho)*norm(yd))) - n*log(norm(yd));
end

function metric = metric_cor(y,m,s)
ys = y(1:m);
metric = real(ys' * s);
end