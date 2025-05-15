%% general CCA stuff
animal = 'ss2';
lag = 0;

load(sprintf('output/%s_%d.mat',animal,lag));

CFA = eval([animal '.CFA']);
RFA = eval([animal '.RFA']);

CFA(isnan(CFA)) = 0;
RFA(isnan(RFA)) = 0;
% Remove zero rows
CFA(all(~CFA,2),:) = [];
RFA(all(~RFA,2),:) = [];
disp(size(CFA,1));

%% do PCA first

[coefCFA,scoreCFA,latentCFA] = pca(CFA');
[coefRFA,scoreRFA,latentRFA] = pca(RFA');

cvar_CFA = cumsum(latentCFA)/sum(latentCFA);
cvar_RFA = cumsum(latentRFA)/sum(latentRFA);

%%
pcadimCFA = find(cvar_CFA>0.95,1);
pcadimRFA = find(cvar_RFA>0.95,1);


[CV1,CV2,r,score1,score2] = canoncorr(scoreCFA(:,1:pcadimCFA),scoreRFA(:,1:pcadimRFA));

%% orthonormal stuff and getting var capture

% but now, orthonormalize CV1 into orthCV1
orth1 = orth(CV1);

Q = scoreCFA(:,1:pcadimCFA)*orth1;
varcap1 = flipud(diag(cov(Q)) / trace(cov(scoreCFA)));
% sum(varcap1);
plot(varcap1)














%%


% 
% 
% 
% Q = scoreCFA(:,1:pcadimCFA)*CV1;    % this is exactly score1
% var(Q(:,1))     % this equals 1
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% % but now, orthonormalize CV1 into orthCV1
% Q = scoreRFA(:,1:pcadimRFA)*GramSchmidt(CV2);     % this is exactly score1
% var(Q(:,1))     % this is not 1 anymore
% varcap2 = diag(cov(Q)) / trace(cov(scoreRFA))
% 
% 
% 
% 
% 
% 






%%
disp('weeeeee!');