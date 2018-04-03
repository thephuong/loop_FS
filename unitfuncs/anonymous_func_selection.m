%% Create anonymous functions for data generation and for metric
function [fmetric,fGenData] = anonymous_func_selection(rule,dataType,SIMUPARAMS)
switch (dataType)
    case SIMUPARAMS.CONST_dataType_UNISPHERE
        fGenData = @(n,r,ntest_perbloc) fGenUniVec_cpx_bloc(n,r,ntest_perbloc,-1);
        switch (rule)
            case SIMUPARAMS.CONST_ML_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_uni_cpx_bloc(circshift(y,-tau),m,n,rho,s);
            case SIMUPARAMS.CONST_ML_RULE_NormalApprox_POWER
                fmetric = @(y,tau,m,n,rho,s) fmetric_uni_cpx_bloc_NormalApprox_power(circshift(y,-tau),m,n,rho,s);
            case SIMUPARAMS.CONST_ML_RULE_NormalApprox_DIST
                fmetric = @(y,tau,m,n,rho,s) fmetric_uni_cpx_bloc_NormalApprox_dist(circshift(y,-tau),m,n,rho,s);
            case SIMUPARAMS.CONST_COR_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_cor_cpx_bloc(circshift(y,-tau),m,s);
            otherwise
                error('Rule(%d) not supported.',rule);
        end
    case SIMUPARAMS.CONST_dataType_GAUSSIAN
        fGenData = @(n,r,ntest_perbloc) r/sqrt(2*n)*(randn(n,ntest_perbloc) + 1i*randn(n,ntest_perbloc));
        switch (rule)
            case SIMUPARAMS.CONST_ML_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_Normal_cpx_bloc(circshift(y,-tau),m,rho,s);
            case SIMUPARAMS.CONST_COR_RULE
                fmetric = @(y,tau,m,n,rho,s) fmetric_cor_cpx_bloc(circshift(y,-tau),m,s);
            otherwise
                error('Rule(%d) not supported.',rule);
        end
    otherwise
        error('dataType(%d) not supported.',dataType);
end
end