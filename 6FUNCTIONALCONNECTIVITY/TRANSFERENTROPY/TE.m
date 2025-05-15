%% load TE data
homepath = 'C:\Users\mirilab\Documents\Mark\TransferEntropy\';
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
%%
P = get_pvals(observedTE,NULLTE,ALLOWED);

%%
plot_sig(P,[nRFA nCFA]);

%% separate the matrix into CFA->RFA and RFA->CFA | this only works if RFA comes before CFA on the matrix
pCtoR = P(nRFA+1:end, 1:nRFA);
pRtoC = P(1:nRFA, nRFA+1:end);
obsCtoR = observedTE(nRFA+1:end, 1:nRFA);
obsRtoC = observedTE(1:nRFA, nRFA+1:end);

sigRtoC = pRtoC<0.05;
sigCtoR = pCtoR<0.05;

qTE_RtoC = obsRtoC(sigRtoC);
qTE_CtoR = obsCtoR(sigCtoR);

fprintf('%d, %d \n', sum(qTE_RtoC), length(qTE_RtoC));
fprintf('%d, %d \n', sum(qTE_CtoR), length(qTE_CtoR));


%%
sig = P<0.05;
sigTE = observedTE(sig);

figure;
hold on
histogram(observedTE(:), 'BinWidth',1e-5);
histogram(sigTE, 'BinWidth',1e-5);
hold off



%%