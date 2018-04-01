function perrp = fpredict_perr_Normal_margin_corr_cpx(s,m,n,rho)
N = m+n;
norms = norm(s);
norms2 = norms^2;

perrp = zeros(1,N-1);

for tau = 1:m-1
    s2_m_tau = s(end-m+tau+1:end);   %last m-tau elements
    s1_m_tau = s(1:m-tau);           %first m-tau elements
    s2_tau = s(end-tau+1:end);       %last tau elements
    mu = norms2 - real(s1_m_tau'*s2_m_tau);
    sigma2 = mu + norm(s2_tau)^2 * rho/2;
    perrp(tau) = normcdf(0,mu,sqrt(sigma2));
end

for tau = -m+1:-1
    atau = abs(tau);
    s1_m_atau = s(1:m-atau);
    s1_atau = s(1:atau);
    s2_m_atau = s(end-m+atau+1:end);
    mu = norms2 - real(s1_m_atau'*s2_m_atau);
    sigma2 = mu + norm(s1_atau)^2 * rho/2;
    perrp(tau) = normcdf(0,mu,sqrt(sigma2));
end

perrp(m:end-m+1) = normcdf(0, norms2, norms * sqrt(1+rho/2));