function pp = fdoubly_ncF_bloc(m1,lambda1,m2,lambda2,z0)
    NTEST = 1e5;
	ntest_perbloc=1e3;
	nbloc = ceil(NTEST/ntest_perbloc);
	nn = 0;
	for ibloc = 1:nbloc
		nn = nn + sum((ncx2rnd(m1,lambda1,ntest_perbloc,1) < z0 * ncx2rnd(m2,lambda2,ntest_perbloc,1)));
		if (nn > 100); break; end
	end
	pp = nn/ibloc/ntest_perbloc;
end