clear all

n = 10;
m = randi(n);

rhodB = 0;
rho = 10^(rhodB/10);
r = sqrt(n*rho);

s = 2*(rand(m,1) > 0.5) - 1;
%s = randn(m,1);
%s = ones(m,1);
% s = [5; zeros(m-1,1)];
s = s / norm(s);

NTEST = 1e6;

histo = zeros(NTEST,1);
histoo = zeros(NTEST,1);
ndm = zeros(NTEST,1);
for i = 1:NTEST
	d = fGenUniVec_cpx(n,r);
	dm = d(1:min(n,m));
	ndm(i) = norm(dm);
	histo(i) = real(s' * dm);
    histoo(i) = real(d(1));
end

figure;
hh = histogram(histo);
title(sprintf('Histogram of D_{0r} n=%d rho=%ddB',n,rhodB));

DISP_FACTOR=2;
xx = hh.BinLimits(1):DISP_FACTOR*hh.BinWidth:hh.BinLimits(2);
pxx = 1/norm(s)/sqrt(n*rho) * 1/beta(n-0.5,0.5) * (1-xx.^2/(norm(s)^2*n*rho)).^(n-1.5);
hold on;
plot(xx(1:2:end),pxx(1:2:end)*hh.BinWidth*NTEST,'ro','LineWidth',1);
%normal approx
pnx = normpdf(xx,0,norm(s)*sqrt(n*rho/(2*n-3)));
plot(xx(2:2:end),pnx(2:2:end)*hh.BinWidth*NTEST,'m+','LineWidth',2);
legend('histogram','Beta PDF','Normal approximation');

% figure;
% histogram(histoo,200);
% title(sprintf('Histogram of 1 coordinate n=%d m=%d rho=%d',n,m,rhodB));

% figure;
% histogram(ndm,200);
% title("Histogram of norm of dm");

% figure;
% plot(fDistCoorUniNSphere(n));