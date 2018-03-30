%% Normal approximation of epsilon^*(k,n,rho)
% k: number of information bit
% n: codeword length (complex)
% rho: LINEAR SNR
function eD = epsilon_D_cpx(k,n_tab,rhoD_tab)
    nreal = 2*n_tab(:);
    rhoD_tab = rhoD_tab(:); % noise is normalized 0.5 thus rho keeps intact
    % capa
    C = 0.5*log2(1+rhoD_tab);
    % dispersion
    V = rhoD_tab .* (2+rhoD_tab)/2 ./ (1+rhoD_tab).^2;
    % epsilon_D
    eD = qfunc((nreal.*C + 0.5*log2(nreal) - k) ./ sqrt(nreal.*V) / log2(exp(1)));
end