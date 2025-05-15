%% pca and pls/cca for RFA 

figure;
subplot(1,2,1);
hold on
errorbar(1:ncomps,cumsum(meanPCA_RFA),semPCA_RFA,'-k')
errorbar([1:ncomps]+0.5,cumsum(meanPLS_RFA),semPLS_RFA,'-m')
ylim([0 1])
title('PCA and PLS')

subplot(1,2,2)
hold on
errorbar(1:ncomps,cumsum(meanPCA_RFA),semPCA_RFA,'-k')
errorbar([1:ncomps]+0.5,cumsum(meanCCA_RFA2),semCCA_RFA2,'-m')
ylim([0 1])
title('PCA and CCA')