%% PLS no lag, compare varcap to PCA

% after running PLSSVD the struct D should have everything
% the num_PCs comes from assuming you've run PLSSVD before which defines it there
ncomps = num_PCs;

PCA_varcapRFA = [];
PCA_varcapCFA = [];
PLS_varcapRFA = [];
PLS_varcapCFA = [];

COVS = [];
for i=1:length(D)
    disp(i);
    lags = D(i).lags;
    lag0 = find(lags==0);
    RFA = D(i).RFA_RG;
    CFA = D(i).CFA_RG(:,:,lag0);
    

    [coefRFA,scoreRFA,latentRFA] = pca(RFA');
    [coefCFA,scoreCFA,latentCFA] = pca(CFA');
    vRFA = latentRFA / sum(latentRFA);
    vCFA = latentCFA / sum(latentCFA);

    vRFA = vRFA(1:ncomps);
    vCFA = vCFA(1:ncomps);

    [~,~,covs,vcRFApls,vcCFApls] = pls_svd(RFA',CFA');
    vcRFApls = vcRFApls(1:ncomps);
    vcCFApls = vcCFApls(1:ncomps);
    covs = covs(1:ncomps);    



    PCA_varcapRFA = [PCA_varcapRFA; vRFA'];
    PCA_varcapCFA = [PCA_varcapCFA; vCFA'];
    PLS_varcapRFA = [PLS_varcapRFA; vcRFApls];
    PLS_varcapCFA = [PLS_varcapCFA; vcCFApls];
    COVS = [COVS; covs'];
end

%%
n = length(D);

meanPCA_RFA = mean(PCA_varcapRFA);
meanPCA_CFA = mean(PCA_varcapCFA);
meanPLS_RFA = mean(PLS_varcapRFA);
meanPLS_CFA = mean(PLS_varcapCFA);

semPCA_RFA = std(PCA_varcapRFA) / sqrt(n);
semPCA_CFA = std(PCA_varcapCFA) / sqrt(n);
semPLS_RFA = std(PLS_varcapRFA) / sqrt(n);
semPLS_CFA = std(PLS_varcapCFA) / sqrt(n);

mean_COVS = mean(COVS);
sem_COVS = std(COVS) / sqrt(n);

%%

figure;
subplot(2,1,1);
hold on
boundedline(1:ncomps,cumsum(meanPCA_RFA),semPCA_RFA,'-k','alpha')
boundedline(1:ncomps,cumsum(meanPLS_RFA),semPLS_RFA,'-m','alpha')
hold off
axis tight
ylim([0 1])
ylabel('cumulative variance captured')


subplot(2,1,2);
hold on
boundedline(1:ncomps,cumsum(meanPCA_CFA),semPCA_CFA,'-k','alpha')
boundedline(1:ncomps,cumsum(meanPLS_CFA),semPLS_CFA,'-c','alpha')
hold off
axis tight
ylim([0 1])
xlabel('component')
ylabel('cumulative variance captured')



%%


figure;
subplot(2,1,1);
hold on
errorbar(1:ncomps,cumsum(meanPCA_RFA),semPCA_RFA,'-k')
errorbar(1:ncomps,cumsum(meanPLS_RFA),semPLS_RFA,'-m')
hold off
axis tight
ylim([0 1])
ylabel('cumulative variance captured')


subplot(2,1,2);
hold on
errorbar(1:ncomps,cumsum(meanPCA_CFA),semPCA_CFA,'-k')
errorbar(1:ncomps,cumsum(meanPLS_CFA),semPLS_CFA,'-c')
hold off
axis tight
ylim([0 1])
xlabel('component')
ylabel('cumulative variance captured')




%% find what percent the PLS components capture of the PCA variance

% first 10 CVs
[meanpctVC_CFA,sempctVC_CFA] = find_pct_varcap(PCA_varcapCFA,PLS_varcapCFA);
[meanpctVC_RFA,sempctVC_RFA] = find_pct_varcap(PCA_varcapRFA,PLS_varcapRFA);









%%