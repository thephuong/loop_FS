function d = fGenUniVec_cpx(n,r)
w = 1/sqrt(2) * (randn(n,1) + 1i * randn(n,1));
d = r/norm(w) * w;