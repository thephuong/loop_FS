function ny = fvecwisenorm(y)
if (version('-release')=='2017b')
    ny = vecnorm(y);
else
    ny = sqrt(sum(y.*conj(y),1));
end