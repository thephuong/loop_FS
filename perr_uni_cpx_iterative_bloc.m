function perr = perr_uni_cpx_iterative_bloc(m,n,s,rho,rule,NTEST,nbiteration,nbframe)

if nargin < 6; NTEST = 1e6; end
if nargin < 5; rule = 2; NTEST = 1e6; end
if nargin < 7; nbiteration = 2; end
N = m+n;
r = sqrt(n*rho);

ntest_perbloc = 1000;
nbloc = ceil(NTEST/ntest_perbloc);
sbloc = repmat(s,1,ntest_perbloc);
% nerr = zeros(NTEST,1);
nerrs = 0;
for ibloc = 1:nbloc
    x = zeros(nbframe*N,ntest_perbloc);
    for iframe = 1:nbframe
        x((iframe-1)*N+1:(iframe-1)*N+N,:) = [sbloc;fGenUniVec_cpx_bloc(n,r,ntest_perbloc)];
    end
    tau0 = randi(nbframe*N) - 1;
%     y = x + 1/sqrt(2)*(randn(nbframe*N,ntest_perbloc)+1i*randn(nbframe*N,ntest_perbloc));
    y = circshift(x + 1/sqrt(2)*(randn(nbframe*N,ntest_perbloc)+1i*randn(nbframe*N,ntest_perbloc)),tau0);
    
%     ptau = (1/N)*ones(N,ntest_perbloc);
%     for iter = 1:nbiteration
%         for iframe = 1:nbframe
%             ptau = fmetric_alltau_bloc(y((iframe-1)*N+1:(iframe-1)*N+N,:),m,n,rho,s,ptau,rule);
%         end
%     end
%     nerrs = nerrs + fis_err_bloc(ptau, mod(tau0,N)+1);
    logptau = -log(N)*ones(N,ntest_perbloc);
    for iter = 1:nbiteration
        for iframe = 1:nbframe
            logptau = fmetric_log_alltau_bloc(y((iframe-1)*N+1:(iframe-1)*N+N,:),m,n,rho,s,logptau,rule);
        end
    end
    nerrs = nerrs + fis_err_bloc(logptau, mod(tau0,N)+1);

    if (nerrs > 500)
        break;
    end
end

% perr = mean(nerr);
perr = nerrs/ibloc/ntest_perbloc;
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

function logpMAP = fmetric_log_alltau_bloc(y,m,n,rho,s,logptau,rule)
N = m+n;
temp = zeros(N,size(logptau,2));
if (rule == 1)
    temp(1,:) = fmetric_uni_cpx_bloc(y,m,n,rho,s);
    for tau=1:N-1
        temp(tau+1,:) = fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s);
    end
else %rule==2
    temp(1,:) = metric_cor2(y,m,s);
    for tau=1:N-1
        temp(tau+1,:) = metric_cor2(circshift(y,-tau),m,s);
    end
end
logpMAP = temp + logptau;
end

function pMAP = fmetric_alltau_bloc(y,m,n,rho,s,ptau,rule)
N = m+n;
temp = zeros(N,size(ptau,2));
if (rule == 1)
    temp(1,:) = fmetric_uni_cpx_bloc(y,m,n,rho,s);
    for tau=1:N-1
        temp(tau+1,:) = fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s);
    end
else %rule==2
    temp(1,:) = metric_cor2(y,m,s);
    for tau=1:N-1
        temp(tau+1,:) = metric_cor2(circshift(y,-tau),m,s);
    end
end
pMAP_unnormalized = temp + log(ptau);
pMAP = exp(pMAP_unnormalized) ./ repmat(sum(exp(pMAP_unnormalized),1),N,1);
% pMAP = temp;
end

function metric = metric_cor2(y,m,s)
metric = 2*real(s'*y(1:m,:));
end