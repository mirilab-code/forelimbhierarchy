%% FA
clear all

load('03312020_binned_spike_counts.mat');


runIdx = 1;           % Results will be saved in baseDir/mat_results/runXXX/, where
% XXX is runIdx. Use a new runIdx for each dataset.
baseDir = '.';        % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;     % Set to true to save train and test data (not recommended)
binWidth = 10;        % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 5;         % Number of cross-validation folds (0 means no cross-validation)
xDims = {0,0};       % The number of latents for each group
yDims = [56 81];      % Number of observed features (neurons) in each group (area)
maxIters = 1e8;       % Limit the number of EM iterations (not recommended, in general)
randomSeed = 0;       % Seed the random number generator, for reproducibility
numGroups = length(yDims); % Number of observation groups in the data

fit_fa(runIdx, all_seq, ...
    'baseDir', baseDir, ...
    'binWidth', binWidth, ...
    'numFolds', numFolds, ...
    'xDims', xDims, ...
    'yDims', yDims, ...
    'maxIters', maxIters, ...
    'parallelize', false, ...  % Only relevant if cross-validating
    'randomSeed', randomSeed, ...
    'overwriteExisting', overwriteExisting, ...
    'saveData', saveData);



%% DLAG

runIdx = 1;               % Results will be saved in baseDir/mat_results/runXXX/,
% where XXX is runIdx. Use a new runIdx for each dataset.
baseDir = '.';            % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;         % Set to true to save train and test data (not recommended)
method = 'dlag';          % For now this is the only option, but that may change in the near future
binWidth = 20;            % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;             % Number of cross-validation folds (0 means no cross-validation)
xDims_across = 4;         % This number of across-group latents matches the synthetic ground truth
xDims_within = {2, 2};    % These numbers match the within-group latents in the synthetic ground truth
yDims = [73 116];          % Number of observed features (neurons) in each group (area)
rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
segLength = 25;           % Largest trial segment length, in no. of time points
init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
learnDelays = true;       % Set to false if you want to fix delays at their initial value
maxIters = 5e3;           % Limit the number of EM iterations (not recommended for final fitting stage)
freqLL = 10;              % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 100;          % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.01;        % Private noise variances will not be allowed to go below this value
verbose = true;           % Toggle printed progress updates
randomSeed = 0;           % Seed the random number generator, for reproducibility

fit_dlag(runIdx, all_seq, ...
    'baseDir', baseDir, ...
    'method', method, ...
    'binWidth', binWidth, ...
    'numFolds', numFolds, ...
    'xDims_across', xDims_across, ...
    'xDims_within', xDims_within, ...
    'yDims', yDims, ...
    'rGroups', rGroups,...
    'startTau', startTau, ...
    'segLength', segLength, ...
    'init_method', init_method, ...
    'learnDelays', learnDelays, ...
    'maxIters', maxIters, ...
    'freqLL', freqLL, ...
    'freqParam', freqParam, ...
    'minVarFrac', minVarFrac, ...
    'parallelize', false, ... % Only relevant for cross-validation
    'verbose', verbose, ...
    'randomSeed', randomSeed, ...
    'overwriteExisting', overwriteExisting, ...
    'saveData', saveData);