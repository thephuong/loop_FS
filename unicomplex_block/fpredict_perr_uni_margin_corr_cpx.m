function perrp = fpredict_perr_uni_margin_corr_cpx(s,m,n,rho)
N = m+n;
norms = norm(s);
norms2 = norms^2;

perrp = zeros(1,N-1);
for tau = 1:m-1
%     s1tau = s(1:tau);               %first tau elements
    s2m_tau = s(end-m+tau+1:end);   %last m-tau elements
    s1m_tau = s(1:m-tau);           %first m-tau elements
    s2tau = s(end-tau+1:end);       %last tau elements
    mu = norms2 - real(s2m_tau' * s1m_tau);
    sigma2 = mu + norm(s2tau)^2 * rho * n/(2*n-3);
    perrp(tau) = normcdf(0,mu,sqrt(sigma2));
end
for tau = N-1-(m-1)+1:N-1
    tauh = tau - n;
    s1tau = s(1:tauh);               %first tau elements
%     s2m_tau = s(end-m+tau+1:end);   %last m-tau elements
    s1m_tau = s(1:m-tauh);           %first m-tau elements
    s2tau = s(end-tauh+1:end);       %last tau elements
    mu = norms2 - real(s2tau' * s1tau);
    sigma2 = mu + norm(s1m_tau)^2 * rho * n/(2*n-3);
    perrp(tau) = normcdf(0,mu,sqrt(sigma2));
end
perrp(m:end-m+1) = normcdf(0, norms^2, norms * sqrt(1+n*rho/(2*n-3)));