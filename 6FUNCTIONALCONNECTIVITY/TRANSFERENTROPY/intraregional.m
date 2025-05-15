%% check intra regional connectivity
datapath = 'data\TEresults10kpairs';
addpath(datapath);

files = fullfile(datapath, '0327*.mat');
matFiles = dir(files);
for i = 1:length(matFiles)
    baseFileName = fullfile(datapath, matFiles(i).name);
    load(baseFileName);
end
disp('done loading!');

NULLTE = cat(3,nullTE,nullTE2);
ALLOWED = cat(3,allowed_shifts,allowed_shifts2);

P = get_pvals(observedTE,NULLTE,ALLOWED);

%%
plot_sig(P,[nRFA nCFA]);

%% this only works if RFA comes before CFA on the matrix
pCtoR = P(nRFA+1:end, 1:nRFA);
pRtoC = P(1:nRFA, nRFA+1:end);
pRtoR = P(1:nRFA,1:nRFA);
pCtoC = P(nRFA+1:end,nRFA+1:end);
obsCtoR = observedTE(nRFA+1:end, 1:nRFA);
obsRtoC = observedTE(1:nRFA, nRFA+1:end);
obsRtoR = observedTE(1:nRFA,1:nRFA);
obsCtoC = observedTE(nRFA+1:end,nRFA+1:end);

%%
alpha = 0.05;
bw = 0.01;

pvalsRtoC = pRtoC(:);
pvalsCtoR = pCtoR(:);
pvalsRtoR = pRtoR(:);
pvalsCtoC = pCtoC(:);


figure;
subplot(2,3,[1,3])
% RFA to CFA
hold on
histogram(pvalsRtoC,'BinWidth',bw);
histogram(pvalsCtoR,'BinWidth',bw);
legend({'RFA to CFA','CFA to RFA'});

subplot(2,3,4);
histogram(pvalsRtoR,'BinWidth',bw);
title('intraregional RFA')

subplot(2,3,5);
histogram(pvalsCtoC,'BinWidth',bw);
title('intraregional CFA')

subplot(2,3,6);
hold on
histogram(pvalsRtoR,'BinWidth',bw, 'Normalization','probability');
histogram(pvalsCtoC,'BinWidth',bw, 'Normalization','probability');
hold off
legend({'intraregional RFA','intraregional CFA'})












%%