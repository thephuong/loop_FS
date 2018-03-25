function metric = fmetric_uni_cpx_bloc(y,m,n,rho,s)
% ys = y(1:m,:);
normyd = fvecwisenorm(y(m+1:end,:));
metric = 2*real(s'*y(1:m,:)) + log(besseli(n-1,2*sqrt(n*rho)*normyd)) - n*log(normyd);
end