%% Generate data D uniformly distributed on a complex sphere
% embbed in n dimensional complex space
% radius r
% generate ntest_perbloc vector D
% tau_frac: fraction factor of two D vectors
function d = fGenUniVec_cpx_bloc(n,r,ntest_perbloc,varargin)
if isempty(varargin)
    tau_frac = 0;
else
    tau_frac = varargin{1};
end
w = 1/sqrt(2) * (randn(n,ntest_perbloc) + 1i * randn(n,ntest_perbloc));
normw = repmat(fvecwisenorm(w),n,1);
d = r * w ./ normw;
switch tau_frac
    case 0
        return;
    case -1
        % Random fraction
        w1 = 1/sqrt(2) * (randn(n,ntest_perbloc) + 1i * randn(n,ntest_perbloc));
        d1 = r * w1 ./ repmat(fvecwisenorm(w1),n,1);
        for ivec = 1:ntest_perbloc
            ttau_frac = randi(n)-1;
            d(n-ttau_frac+1:n,ivec) = d1(n-ttau_frac+1:n,ivec);
        end
    otherwise
        w1 = 1/sqrt(2) * (randn(n,ntest_perbloc) + 1i * randn(n,ntest_perbloc));
        d1 = r * w1 ./ repmat(fvecwisenorm(w1),n,1);
        d = [d(1:n-tau_frac,:) ; d1(n-tau_frac+1:n,:)];
end
% d = zeros(n,ntest_perbloc);
% for i = 1:ntest_perbloc
%     w = 1/sqrt(2) * (randn(n,1) + 1i * randn(n,1));
%     d(:,i) = r/norm(w) * w;
% end