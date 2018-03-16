function d = fGenUniVec_cpx_bloc(n,r,ntest_perbloc)
w = 1/sqrt(2) * (randn(n,ntest_perbloc) + 1i * randn(n,ntest_perbloc));
normw = repmat(fvecwisenorm(w),n,1);
d = r * w ./ normw;
% d = zeros(n,ntest_perbloc);
% for i = 1:ntest_perbloc
%     w = 1/sqrt(2) * (randn(n,1) + 1i * randn(n,1));
%     d(:,i) = r/norm(w) * w;
% end