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
parpool('local',48);
dataSet = 3;
% load spike trains
switch dataSet
	% simple network consisting of two units with one connection and
	% non-stationarity in their spike rate
	case 1
		D = load('../../Results/ExampleDataSet_2Units_Exogenous.mat');
		Topology = D.Topology;
        time_Topology_ms = D.time_Topology_ms;
        SpikeTrains = D.SpikeTrains;
        sample_Hz = D.sample_Hz;
        sampleTime_ms = 1000 / sample_Hz;
		
	% monkey data - Fig. 7 and S1 of the paper by Casile et al.
	case 2
		D = load('../../Results/ExampleDataSet_MonkeyData.mat');
        SpikeTrains = D.SpikeTrains;
        sample_Hz = D.sample_Hz;
        sampleTime_ms = 1000 / sample_Hz;

    case 3
        D = load('/projects/p31350/Granger_Quest/10242020_granger40_newWin.mat');
        %SpikeTrains = D.high_FR_granger([7,18,2,40,43,56],750:1050,:);
        %SpikeTrains = D.high_FR_granger(1:6,750:1050,:);
        SpikeTrains = D.high_FR_granger_comb(:,:,:);
        SpikeTrains(:,end,:) = [];
        sample_Hz = 1000;
        sampleTime_ms = 1000 / sample_Hz;
end
%%
%clear D

a = zeros(1,size(SpikeTrains,1));
for ii = 1:size(SpikeTrains,1)
    a(ii) = length(find(SpikeTrains(ii,:,:)==1));
end
disp(a)
% SpikeTrains= SpikeTrains(find(a>20),:,:);
% b = zeros(1,size(SpikeTrains,1));
% for ii = 1:size(SpikeTrains,1)
%     b(ii) = length(find(SpikeTrains(ii,:,:)==1));
% end
% disp(b')

%%

% get information about spike trains
[nNeurons, lenTrial_samples, nTrials] = size(SpikeTrains);
lenTrial_ms = lenTrial_samples * sampleTime_ms;
disp('Done with the spikes trains! Now I fit GLM models');

% ---------------- define global regressor ------------------------
% number of windows used to divide each trial
globalRegressor.nBins = [1, 5, 10, 15, 20, 25, 30];
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
historyRegressorNBins = [2:2:historyRegressor.maxNBins];
%%
% now run the Granger causality method
OutStruct = runGranger_G_ETM(SpikeTrains, globalRegressor, historyRegressor, historyRegressorNBins);

% now save all the results
fName = ['./10242020_newWin_Out.mat'];
fprintf('Saving %s \n', fName);
save(fName);
%%
% ... and now let's plot the result


