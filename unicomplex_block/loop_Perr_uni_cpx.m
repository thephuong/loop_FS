clear all; close all;
alpha_tab = -2:1;

for ii = 1:length(alpha_tab)
    alpha = alpha_tab(ii);
    fprintf('Alpha=%d\n',alpha);
    test_Perr_uni_cpx_alpha;
    fprintf('Alpha=%d DONE ...\n',alpha);
%     close;close;
%     close all;
end