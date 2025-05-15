
% after running PLSSVD the struct D should have everything
% the num_PCs comes from assuming you've run PLSSVD before which defines it there
ncomps = num_PCs;

PCA_varcapRFA = [];
PCA_varcapCFA = [];
CCA_varcapRFA = [];
CCA_varcapCFA = [];

CORRS = [];
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

    [~,~,r,RFAcca,CFAcca,vcRFA,vcCFA] = pca_then_cca(RFA,CFA,ncomps);
  



    PCA_varcapRFA = [PCA_varcapRFA; vcRFA.pca'];
    PCA_varcapCFA = [PCA_varcapCFA; vcCFA.pca'];
    CCA_varcapRFA = [CCA_varcapRFA; vcRFA.cca'];
    CCA_varcapCFA = [CCA_varcapCFA; vcCFA.cca'];
    CORRS = [CORRS; r];
end

%%
n = length(D);

meanPCA_RFA2 = mean(PCA_varcapRFA);
meanPCA_CFA2 = mean(PCA_varcapCFA);
meanCCA_RFA2 = mean(CCA_varcapRFA);
meanCCA_CFA2 = mean(CCA_varcapCFA);

semPCA_RFA2 = std(PCA_varcapRFA) / sqrt(n);
semPCA_CFA2 = std(PCA_varcapCFA) / sqrt(n);
semCCA_RFA2 = std(CCA_varcapRFA) / sqrt(n);
semCCA_CFA2 = std(CCA_varcapCFA) / sqrt(n);

mean_CORRS2 = mean(CORRS);
sem_CORRS2 = std(CORRS) / sqrt(n);

%%

figure;
subplot(2,1,1);
hold on
boundedline(1:ncomps,cumsum(meanPCA_RFA2),semPCA_RFA2,'-k','alpha')
boundedline(1:ncomps,cumsum(meanCCA_RFA2),semCCA_RFA2,'-m','alpha')
hold off
axis tight
ylim([0 1])
ylabel('cumulative variance captured')


subplot(2,1,2);
hold on
boundedline(1:ncomps,cumsum(meanPCA_CFA2),semPCA_CFA2,'-k','alpha')
boundedline(1:ncomps,cumsum(meanCCA_CFA2),semCCA_CFA2,'-c','alpha')
hold off
axis tight
ylim([0 1])
xlabel('component')
ylabel('cumulative variance captured')



%%


figure;
subplot(3,1,1);
hold on
errorbar(1:ncomps,cumsum(meanPCA_RFA2),semPCA_RFA2,'-k')
errorbar(1:ncomps,cumsum(meanCCA_RFA2),semCCA_RFA2,'-m')
hold off
axis tight
ylim([0 1])
ylabel('cumulative variance captured')


subplot(3,1,2);
hold on
errorbar(1:ncomps,cumsum(meanPCA_CFA2),semPCA_CFA2,'-k')
errorbar(1:ncomps,cumsum(meanCCA_CFA2),semCCA_CFA2,'-c')
hold off
axis tight
ylim([0 1])
ylabel('cumulative variance captured')

subplot(3,1,3);
errorbar(mean_CORRS2,sem_CORRS2,'-k')
axis tight
ylim([0 1])
xlabel('component')


%% find what percent the CVs capture of the PCA variance
% first 10 CVs
[meanpctVC_CFA,sempctVC_CFA] = find_pct_varcap(PCA_varcapCFA,CCA_varcapCFA);
[meanpctVC_RFA,sempctVC_RFA] = find_pct_varcap(PCA_varcapRFA,CCA_varcapRFA);



%%





%% first 10 CVs%
% cs_cca_CFA = cumsum(meanCCA_CFA2);
% cs_pca_CFA = cumsum(meanPCA_CFA2);
% pct_varcap_CFA = cs_cca_CFA ./ cs_pca_CFA;
% pct_varcap_CFA10 = pct_varcap_CFA(10)
% 
% cs_cca_RFA = cumsum(meanCCA_RFA2);
% cs_pca_RFA = cumsum(meanPCA_RFA2);
% pct_varcap_RFA = cs_cca_RFA ./ cs_pca_RFA;
% pct_varcap_RFA10 = pct_varcap_RFA(10)





%%
