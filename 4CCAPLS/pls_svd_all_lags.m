function [sum_cov]  = pls_svd_all_lags(M1,M2_lags,lags,ncomps)


sum_cov = [];
for i=1:length(lags)
    X = M1;
    Y = M2_lags(:,:,i);

    [~,~,covs,~,~] = pls_svd(X',Y');
    s = sum(covs(1:ncomps));
    sum_cov = [sum_cov s];


end






end