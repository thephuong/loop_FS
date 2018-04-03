function metric = fmetric_uni_cpx_bloc_NormalApprox_dist(y,m,n,rho,s)
rho_Normal = rho*2*n/(2*n-3);
ys = y(1:m,:);
metric = real(s'*ys) - rho_Normal/2/(rho_Normal+1)*(fvecwisenorm(ys)).^2;
end