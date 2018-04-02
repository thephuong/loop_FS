function [s,rho] = fGenSyncWord(m,n,rho_tot,alpha,stype)
    N = m+n;
    rho = N/(m*alpha+n)*rho_tot; %rho_d
    rho_s = alpha*rho; %rho_s

    if (stype == 1)
        s = sqrt(rho_s)*fGenUniVec_cpx(m,sqrt(m*rho));
    elseif (stype == 2)
        %Zadoff Chu, force m impair
        s = sqrt(rho_s)*lteZadoffChuSeq(1,m);
    else
        s = sqrt(rho_s)*(2*(rand(m,1) > 0.5)-1);
    end
end