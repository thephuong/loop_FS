function perrp = fpredict_perr_uni_margin_ML_cpx(s,m,n,rho,ESTIM_SMALL_TAU)

if nargin < 5; ESTIM_SMALL_TAU=0; end

N = m+n;
norms = norm(s);
norms2 = norms^2;

%pdf of z = 2*norm(y)^2
lambda = 2*(n*rho+norms2);
mu = 2*N + lambda;
stdist = sqrt(2*(2*N+2*lambda));
dz = 1;
zz = max(floor(mu-3*stdist),0):dz:ceil(mu+3*stdist);
pzz = ncx2pdf(zz,2*N,lambda);
Petau_yy = zeros(size(zz));
sigma2 = 1/2+rho*n/(2*n-3);
parfor iy2 = 1:length(zz)
    yy = zz(iy2)/2;
    A = sqrt(n*rho/yy) - (2*n+1)/4/yy;
    Petau_yy(iy2) = fdoubly_ncF_bloc(2*m,norms2/A^2/sigma2,2*m,2*(1-1/A)^2*norms2,0.5/sigma2);
end
perrp = ones(1,N-1) * sum(Petau_yy .* (pzz*dz));

if (ESTIM_SMALL_TAU == 0)
    return
end

% For 0<tau<m
alpha3 = 1/2+rho*n/(2*n-3);
for tau = 1:m-1
    s1_atau = s(1:tau);
    s1_m_atau = s(1:m-tau);
    s2_m_atau = s(tau+1:m);
    s2_atau = s(m-tau+1:m);
    Petau_yy = zeros(size(zz));
    parfor iy2 = 1:length(zz)
        yy = zz(iy2)/2;
        A = sqrt(n*rho/yy) - (2*n+1)/4/yy;
        lambda1 = 2*(1-1/A)^2*norm(s1_atau)^2;
        lambda3 = norm(s2_atau)^2/A^2/alpha3;
        mu2 = (norm(s1_m_atau)^2 - norm(s2_m_atau)^2)/A^2 + 2/A * real((s2_m_atau - s1_m_atau)'*s2_m_atau);
        sqrtvar2 = norm(s2_m_atau-s1_m_atau)*sqrt(2)/A;
        Petau_yy(iy2) = linear_combin_tau_lt_m(1/2,lambda1,mu2,sqrtvar2,alpha3,lambda3,tau);
    end
    perrp(tau) = sum(Petau_yy .* (pzz*dz));
    if (ESTIM_SMALL_TAU == 2)
        fprintf('tau=%d perr=%.3e\n',tau,perrp(tau));
    end
end

for tau = -m+1:-1
    atau = abs(tau);
    s1_atau = s(1:atau);
    s1_m_atau = s(1:m-atau);
    s2_m_atau = s(atau+1:m);
    s2_atau = s(m-atau+1:m);
    Petau_yy = zeros(size(zz));
    parfor iy2 = 1:length(zz)
        yy = zz(iy2)/2;
        A = sqrt(n*rho/yy) - (2*n+1)/4/yy;
        lambda1 = 2*(1-1/A)^2 * norm(s2_atau)^2;
        lambda3 = norm(s1_atau)^2/A^2/alpha3;
        mu2 = (norm(s2_m_atau)^2 - norm(s1_m_atau)^2)/A^2 + 2/A * real((s1_m_atau - s2_m_atau)'*s1_m_atau);
        sqrtvar2 = norm(s1_m_atau-s2_m_atau)*sqrt(2)/A;
        Petau_yy(iy2) = linear_combin_tau_lt_m(1/2,lambda1,mu2,sqrtvar2,alpha3,lambda3,atau);
    end
    perrp(tau+N) = sum(Petau_yy .* (pzz*dz));
    if (ESTIM_SMALL_TAU == 2)
        fprintf('tau=%d perr=%.3e\n',tau,perrp(tau+N));
    end
end
end

% function pp = doubly_ncF(m1,lambda1,m2,lambda2,z0)
    % NTEST = 1e5;
    % nn = 0;
    % for i = 1:NTEST
        % nn = nn + (ncx2rnd(m1,lambda1) < z0 * ncx2rnd(m2,lambda2));
        % if (nn > 50); break; end
    % end
    % pp = nn/i;
% end

function pp = linear_combin_tau_lt_m(alpha1,lambda1,mu2,sqrtvar2,alpha3,lambda3,tau)
    NTEST = 1e5;
    nn = 0;
    dof = 2*tau;
    for i = 1:NTEST
        nn = nn + (alpha3*ncx2rnd(dof,lambda3)+normrnd(mu2,sqrtvar2)-alpha1*ncx2rnd(dof,lambda1) < 0);
        if (nn > 50); break; end
    end
    pp = nn/i;
end