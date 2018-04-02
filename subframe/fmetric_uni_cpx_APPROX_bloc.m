function metric = fmetric_uni_cpx_APPROX_bloc(y,m,n,rho,s)
normy = fvecwisenorm(y);
A = sqrt(n*rho)./normy - (2*n+1)/4./(normy.^2);
metric = -fvecwisenorm(y(1:m,:) - repmat(1./A,m,1) .* repmat(s,1,size(y,2)));
end