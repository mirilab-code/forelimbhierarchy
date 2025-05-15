%
% This m-file runs the G-ETM method as described in:
%
% Robust point-process Granger causality analysis in presence of exogenous
% temporal modulations and trial-by-trial variability in spike trains.
%
% by Casile A., Faghih R. T. & Brown E. N.
%
%
% Code tested in Matlab R2019B
%
% author:   Antonino Casile
% toninocasile@gmail.com
%

% data set to load
%parpool('local',48);
%addpath('/projects/p31350/Granger_Quest')
% load spike trains

D = load('C://Users/mirilab/Documents/Adam/reaching/10k_granger_trains/04022020/04022020_granger_train.mat');
C = load('C://Users/mirilab/Documents/Adam/reaching/10k_granger_trains/04022020/04022020_final_list.mat')
pairs_used = C.final_list;
SpikeTrains = D.granger_train(:,:,:);
SpikeTrains(:,end,:) = [];
sample_Hz = 1000;
sampleTime_ms = 1000 / sample_Hz;

a = zeros(1,size(SpikeTrains,1));
for ii = 1:size(SpikeTrains,1)
    a(ii) = length(find(SpikeTrains(ii,:,:)==1));
end
disp(a)

%%

% get information about spike trains
[nNeurons, lenTrial_samples, nTrials] = size(SpikeTrains);
lenTrial_ms = lenTrial_samples * sampleTime_ms;
disp('Done with the spikes trains! Now I fit GLM models');

% ---------------- define global regressor ------------------------
% number of windows used to divide each trial
globalRegressor.nBins = [1, 5, 8, 10, 20, 25];
%globalRegressor.nBins = [10];
% Uncomment the following line and comment out the line above
% to run Kim et al. method
%globalRegressor.nBins = [1];
%%
% temporal duration of each bin of the global regressor
% in BOTH MILLISECONDS and SAMPLES
globalRegressor.binDuration_samples = round(lenTrial_samples ./ globalRegressor.nBins);
globalRegressor.binDuration_ms = globalRegressor.binDuration_samples * sampleTime_ms;

% ---------------- define history regressor -----------------------
% here is the duration of each bin used for the history of the neuron
historyRegressor.binDuration_samples = 3;
historyRegressor.binDuration_ms = historyRegressor.binDuration_samples * sampleTime_ms;

% define maximum number of bins
historyRegressor.maxNBins = 10;
historyRegressor.winHistory_samples = ones(1, historyRegressor.binDuration_samples);
historyRegressor.winHistory_ms = ones(1, historyRegressor.binDuration_ms);

% number of steps for the history regressor that we test
historyRegressorNBins = [10];
%%
% now run the Granger causality method
OutStruct = runGranger_G_ETM_only_half1(SpikeTrains, globalRegressor, historyRegressor, historyRegressorNBins,pairs_used);

% now save all the results
cd('C://Users/mirilab/Documents/Adam/reaching/10k_granger_trains/04022020')
fName = ['./04022020_Out_fixed_only_half1.mat'];
fprintf('Saving %s \n', fName);
save(fName);
%%
% ... and now let's plot the result


