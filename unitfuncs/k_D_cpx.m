%% Normal approximation of log2(M^*(epsilon,n,rho))
% epsilon: max channel decoding error prob.
% n: codeword length (complex)
% rho: LINEAR SNR
function kD = k_D_cpx(epsilon,n_tab,rhoD_tab)
    nreal = 2*n_tab(:);
    rhoD_tab = rhoD_tab(:); % noise is normalized 0.5 thus rho keeps intact
    %capa
    C = 0.5*log2(1+rhoD_tab);
    %dispersion
    V = rhoD_tab .* (2+rhoD_tab)/2 ./ (1+rhoD_tab).^2;
    kD = max(0, nreal .* C - ...
        sqrt(V.*nreal)*log2(exp(1)).*qfuncinv(epsilon) + ...
        0.5*log2(nreal));
end