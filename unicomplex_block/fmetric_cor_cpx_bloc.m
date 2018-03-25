function metric = fmetric_cor_cpx_bloc(y,m,s)
metric = real(s'*y(1:m,:));
end