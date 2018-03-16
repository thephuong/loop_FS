function d = fGen_lteQAM_Vec_cpx_bloc(k,M,n,r,ntest_perbloc)

k = max(k,40);
k = min(k,6144);
switch (M)
    case 1
        modstr = 'BPSK';
    case 2
        modstr = 'QPSK';
    case 4
        modstr = '16QAM';
    case 6
        modstr = '64QAM';
    case 8
        modstr = '256QAM';
    otherwise
        modstr = 'QPSK';
end
sqrtrho = r/sqrt(n);
ratematcher_lenout = n*M;

d = zeros(n,ntest_perbloc);
for itest = 1:ntest_perbloc
    turbo_out = lteTurboEncode(rand(k,1)>0.5);
    ratematcher_out = lteRateMatchTurbo(turbo_out,ratematcher_lenout,0);
    d(:,itest) = sqrtrho * lteSymbolModulate(ratematcher_out, modstr);
end
