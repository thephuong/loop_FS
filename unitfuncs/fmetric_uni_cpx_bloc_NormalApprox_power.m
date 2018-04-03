function metric = fmetric_uni_cpx_bloc_NormalApprox_power(y,m,n,rho,s)
normyd = fvecwisenorm(y(m+1:end,:));
metric = 2*real(s'*y(1:m,:)) + log(besseli(n-1,2*sqrt(n*rho*2*n/(2*n-3))*normyd)) - n*log(normyd);
end