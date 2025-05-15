function [r,pcadimCFA,pcadimRFA,varcapCFA,varcapRFA] = get_CCA_lags(animal,lag)


    load(sprintf('output/%s_%d.mat',animal,lag));

    CFA = eval([animal '.CFA']);
    RFA = eval([animal '.RFA']);

    CFA(isnan(CFA)) = 0;
    RFA(isnan(RFA)) = 0;
    % Remove zero rows
    CFA(all(~CFA,2),:) = [];
    RFA(all(~RFA,2),:) = [];
%     disp(size(CFA,1));

    %% do PCA first

    [coefCFA,scoreCFA,latentCFA] = pca(CFA');
    [coefRFA,scoreRFA,latentRFA] = pca(RFA');

    cvar_CFA = cumsum(latentCFA)./sum(latentCFA);
    cvar_RFA = cumsum(latentRFA)./sum(latentRFA);

    %%
    pcadimCFA = find(cvar_CFA>0.95,1);
    pcadimRFA = find(cvar_RFA>0.95,1);
%     [CV1,CV2,r,score1,score2] = canoncorr(scoreCFA(:,1:pcadimCFA),scoreRFA(:,1:pcadimRFA));
    pcadim = min(pcadimCFA,pcadimRFA);
    [CV1,CV2,r,score1,score2] = canoncorr(scoreCFA(:,1:pcadim),scoreRFA(:,1:pcadim));
    pcadimCFA = pcadim;
    pcadimRFA = pcadim;
    
    %% orthonormal stuff and getting var capture
    orth1 = orth(CV1);
    Q = scoreCFA(:,1:pcadimCFA)*orth1;
    varcapCFA = flipud(diag(cov(Q)) / trace(cov(scoreCFA)));

    orth2 = orth(CV2);
    Q = scoreRFA(:,1:pcadimRFA)*orth2;
    varcapRFA = flipud(diag(cov(Q)) / trace(cov(scoreRFA)));


end










