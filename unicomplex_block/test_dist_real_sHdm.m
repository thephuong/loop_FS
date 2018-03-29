clear all

n = 20;
m = randi(n);

delt = randi(m);
m1 = m-delt;
m2 = delt;

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
% 	d = fGenUniVec_cpx(n,r);
%   dm = d(1:min(n,m));
    d1 = fGenUniVec_cpx(n,r);
    d2 = fGenUniVec_cpx(n,r);
    dm = [d1(n-m2+1:n);d2(1:m1)];

	ndm(i) = norm(dm);
	histo(i) = real(s' * dm);
%     histoo(i) = real(d(1));
end

DISP_FACTOR=1; HISTOBIN=40;
figure;
hh = histogram(histo,HISTOBIN);
title(sprintf('Histogram of D_{0r}, n=%d,m1=%d,m2=%d, rho=%ddB',n,m1,m2,rhodB));
xx = hh.BinLimits(1):DISP_FACTOR*hh.BinWidth:hh.BinLimits(2);
pxx = 1/norm(s)/sqrt(n*rho) * 1/beta(n-0.5,0.5) * (1-xx.^2/(norm(s)^2*n*rho)).^(n-1.5);
hold on;
plot(xx(1:2:end),pxx(1:2:end)*hh.BinWidth*NTEST,'ro','LineWidth',1);
%normal approx
pnx = normpdf(xx,0,norm(s)*sqrt(n*rho/(2*n-3)));
plot(xx(2:2:end),pnx(2:2:end)*hh.BinWidth*NTEST,'m+','LineWidth',2);
legend('histogram','Beta PDF','Normal approximation');

figure;
hh = histogram(histo,HISTOBIN,'Normalization','pdf');
xx = hh.BinLimits(1):DISP_FACTOR*hh.BinWidth:hh.BinLimits(2);
pxx = 1/norm(s)/sqrt(n*rho) * 1/beta(n-0.5,0.5) * (1-xx.^2/(norm(s)^2*n*rho)).^(n-1.5);
pnx = normpdf(xx,0,norm(s)*sqrt(n*rho/(2*n-3)));
title(sprintf('Histogram of D_{0r}, n=%d,m1=%d,m2=%d, rho=%ddB',n,m1,m2,rhodB));
hold on;
plot(xx(1:2:end),pxx(1:2:end),'ro','LineWidth',1);
plot(xx(2:2:end),pnx(2:2:end),'m+','LineWidth',2);
legend('histogram','Beta PDF','Normal approximation');

% figure;
% histogram(histoo,200);
% title(sprintf('Histogram of 1 coordinate n=%d m=%d rho=%d',n,m,rhodB));

% figure;
% histogram(ndm,200);
% title("Histogram of norm of dm");

% figure;
% plot(fDistCoorUniNSphere(n));