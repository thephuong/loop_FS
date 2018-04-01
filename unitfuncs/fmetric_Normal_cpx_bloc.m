function metric = fmetric_Normal_cpx_bloc(y,m,rho,s)
ys = y(1:m,:);
normys = fvecwisenorm(ys);
metric = real(s'*ys) - rho/2/(1+rho) * normys.^2;
end