%% PLS SVD


addpath(genpath('Z:code/'));
%%

matFiles = dir('output/*RFA_CONTROL.mat');  %do this for control
% matFiles = dir('output/*RFA.mat');  %do this for re al calculation
files = {matFiles.name};

for i = 1:length(files)
    fname = files{i};
    d = load(['output/' files{i}]);
    d.session = fname;
    
    nRFA = size(d.RFA_RG,1);
    nCFA = size(d.CFA_RG,1);

%     if(nRFA<20 || nCFA<20)
%         fprintf('skipping %s because of low cell yield \n',fname)
%         continue;
%     end


    fprintf('%s, %d in RFA, %d in CFA \n',fname,nRFA,nCFA);

    if(i==1)
        D = d;
    else
        D = [D d];
    end
end

D = aggregate_animals(D);

lags = D(1).lags;
num_PCs = 25;

%%


%
lag0 = find(lags==0);

cov_lags = [];

for i=1:length(D)
    disp(i);
    RFA = D(i).RFA_RG;
    CFA_lags = D(i).CFA_RG;
    lags = D(i).lags;

    c = pls_svd_all_lags(RFA,CFA_lags,lags,num_PCs);

    cov_lags = [cov_lags; c];
end

%%

mean_sumcovs = mean(cov_lags);
sem_sumcovs = std(cov_lags)/sqrt(size(cov_lags,1));



figure;
subplot(2,1,1);
boundedline(lags,mean_sumcovs,sem_sumcovs, '-b','alpha')
ylabel(sprintf('sum of first %d covariances',num_PCs));

subplot(2,1,2);
boundedline(lags,mean_sumcovs,sem_sumcovs, '-b','alpha')
xlim([-30 30])
ylim([0 700])
ylabel(sprintf('sum of first %d covariances',num_PCs));
xlabel('CFA lag (ms)')






















%%