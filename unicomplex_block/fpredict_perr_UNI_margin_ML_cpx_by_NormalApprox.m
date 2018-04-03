function perrp = fpredict_perr_UNI_margin_ML_cpx_by_NormalApprox(s,m,n,rho)

N = m+n;
norms2 = norm(s)^2;
rho_NormalApprox = rho*2*n/(2*n-3);

lambda1 = 2*(1+rho_NormalApprox)/rho_NormalApprox^2 * norms2;
lambda2 = 2/rho_NormalApprox^2 * norms2;

Pe_tau_gt_m = fdoubly_ncF_bloc(2*m,lambda1, ...
    2*m,lambda2, ...
    1/(1+rho_NormalApprox), 1e9);
perrp = ones(1,N-1) * Pe_tau_gt_m;