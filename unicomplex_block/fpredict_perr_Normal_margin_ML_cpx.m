function perrp = fpredict_perr_Normal_margin_ML_cpx(s,m,n,rho)

N = m+n;
norms2 = norm(s)^2;

lambda1 = 2*(1+rho)/rho^2 * norms2;
lambda2 = 2/rho^2 * norms2;

Pe_tau_gt_m = fdoubly_ncF_bloc(2*m,lambda1, ...
    2*m,lambda2, ...
    1/(1+rho), 1e7);
perrp = ones(1,N-1) * Pe_tau_gt_m;