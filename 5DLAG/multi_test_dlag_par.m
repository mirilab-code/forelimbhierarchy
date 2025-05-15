% =========
% =========
% DLAG DEMO
% =========
%
% This demo shows how we can extract latent variables from multi-population
% data with DLAG (Gokcen et al., 2021). It's recommended to run this script
% section-by-section, rather than all at once (or put a break point before
% Sections 2 and 3, as they may take a long time, depending on your use
% of parallelization).
%
% Section 1 demonstrates how DLAG can be used for exploratory data
% analysis.
%
%     Section 1a fits a DLAG model with a specified number of within-
%     and across-group latent dimensions. Optional arguments are
%     explicitly specified for the sake of demonstration.
%
%     Section 1b takes this model and explores the latent GP timescales and
%     delays. It performs basic inference of within- and across-group
%     latent trajectories. One can compare estimated parameters and
%     trajectories to the ground truth that underlies the demo synthetic
%     data.
%
%     Section 1c demonstrates how to visually scale and order latents
%     according to various metrics, like variance explained within a group
%     or across-group correlation.
%
%     Section 1d demonstrates how to project latents onto ordered sets of
%     modes that capture shared variance explained within a group or
%     correlation across groups.
%
%     Section 1e demonstrates how to denoise observations using a DLAG
%     model. One can compare raw observations to their denoised
%     counterparts.
%
% Section 2 shows how to select optimal DLAG dimensionalities using
% a streamlined cross-validation approach.
%
%     Section 2a estimates the total dimensionality of each group
%     (within + across) by applying FA to each group independently.
%     See example_fa.m for a detailed demo of FA on multi-population data.
%
%     Section 2b determines the optimal DLAG across- and within-group
%     dimensionalities using a streamlined cross-validation approach.
%     The search space is constrained to models such that the across-
%     and within-group dimensionalities for each group add up to the
%     totals established by FA in Section 2b.
%
%     Section 2c fully trains the optimal model selected in Section 2b,
%     assuming the number of EM iterations was limited during
%     cross-validation.
%
% Section 3 demonstrates post-selection inference procedures.
% After selecting the optimal model via cross-validation, these procedures
% can elucidate the uncertainty in parameter estimates.
%
%     Section 3a evaluates how significantly each across-group set of
%     delays deviates from 0, using bootstrapped samples.
%
%     Section 3b constructs bootstrapped confidence intervals for latent
%     delays and timescales. This section involves re-fitting DLAG models
%     to bootstrapped samples, so its runtime is similar to that of
%     Section 2, depending on how much parallelization is used.
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu
%
% Last Revised:
%     26 Feb 2022

%% ================
% 0a) Load demo data
% ===================

cd('Z:\Sarah\DLAG');
addpath(genpath('DLAG-1.0.0'));

clear all
close all

session_path = uigetdir('', 'Data Path');
cd(session_path)
date = session_path(end-7:end);

data = uigetfile('','PCs');
load(data);
work = data(1:end-4);
%% =======================
% 0b) Set up parallelization
% ===========================

% If parallelize is true, all cross-validation folds and bootstrap samples
% will be analyzed in parallel using Matlab's parfor construct.
% If you have access to multiple cores, this provides significant speedup.
parallelize = true;
numWorkers = 8;      % Adjust this to your computer's specs

%% =====================
% 1a) Fitting a DLAG model
% ========================


baseDir = '.';            % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;         % Set to true to save train and test data (not recommended)
method = 'dlag';          % For now this is the only option, but that may change in the near future
binWidth = 20;            % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;             % Number of cross-validation folds (0 means no cross-validation)
yDims = [numCFA numRFA];          % Number of observed features (neurons) in each group (area)
rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
segLength = 20;           % Largest trial segment length, in no. of time points
init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
learnDelays = true;       % Set to false if you want to fix delays at their initial value
maxIters = 5e4;           % Limit the number of EM iterations (not recommended for final fitting stage)
freqLL = 10;              % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 100;          % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.01;        % Private noise variances will not be allowed to go below this value
verbose = false;           % Toggle printed progress updates
randomSeed = 0;           % Seed the random number generator, for reproducibility

%% For loop starts
% parfor (i = 2:6, numWorkers) % Number of across-area latents
%     for x = 4  % Number of within-area latents per area

parfor (x = 2:6, numWorkers)
    for i = 4
        % Let's explicitly define all of the optional arguments, for
        % the sake of demonstration:
        runIdx = x;               % Results will be saved in baseDir/mat_results/runXXX/,
        % where XXX is runIdx. Use a new runIdx for each dataset.
        xDims_across = i;         % This number of across-group latents matches the synthetic ground truth
        xDims_within = {x, x};    % These numbers match the within-group latents in the synthetic ground truth


        fit_dlag(runIdx, DLAG_struct, ...
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

    end
end
%% Save because that part takes forever
cd(session_path)
save(sprintf('%s_new_dlag_model_vw.mat', work)) % Save workspace

% SET BREAKPOINT AFTER THIS



%% =========================================================
% 1b) Explore estimated GP parameters and compare to ground truth
% ================================================================

% 1) Load workspace (va = dim varied; vw = run varied)
% 2) Manually load mat_results file

% Retrieve the fitted model of interest

cd(session_path)
parallelize = true;
numWorkers = 8;

runIdx = xDim_within(1);
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
    'baseDir', baseDir);

% Plot training progress of various quantities. These plots can help with
% troubleshooting, if necessary.
plotFittingProgress(res, ...
    'freqLL', freqLL, ...
    'freqParam', freqParam, ...
    'units', 'ms');

% Plot estimated latents
[seqEst, ~] = exactInferenceWithLL_dlag(DLAG_struct, res.estParams);
plotDimsVsTime_dlag(seqEst, 'xsm', res.estParams, res.binWidth, ...
    'nPlotMax', 1, ...
    'plotSingle', true, ...
    'plotMean', true, ...
    'units', []);

%% ====================================================
% 1c) Visually scale latent trajectories by various metrics
% ==========================================================

% Scale by variance explained
total = false; % true: denominator is total variance; else shared variance
[varexp, ~] = computeVarExp_dlag(res.estParams, total);
[seqEst, sortParams] = scaleByVarExp(seqEst, ...
    res.estParams, ...
    varexp.indiv, ...
    'sortDims', true);
plotDimsVsTime_dlag(seqEst, 'xve', sortParams, res.binWidth, ...
    'nPlotMax', 10, ...
    'plotSingle', true, ...
    'plotMean', false, ...
    'units', []);


%% ==========================================================
% 3a) Evaluate how significantly each set of across-group delays
%     deviates from zero.
%  ================================================================

% Retrieve the best DLAG model
% xDim_across = bestModel.xDim_across;
% xDim_within = bestModel.xDim_within;

res = getModel_dlag(runIdx, xDim_across, xDim_within, 'baseDir', baseDir);

% Save all bootstrap results to a file
boot_fname = generate_inference_fname_dlag(runIdx, ...
    'bootstrapResults', ...
    xDim_across, ...
    xDim_within, ...
    'baseDir',baseDir);
numBootstrap = 1000; % Number of bootstrap samples (the more the better)
% Proportion of bootstrap samples in which the zero-delay model performs
% at least as well as the original model
[delaySig, sigDists] = bootstrapDelaySignificance(DLAG_struct, ...
    res.estParams, ...
    numBootstrap, ...
    'parallelize', parallelize, ...
    'numWorkers', numWorkers);
% Label each delay as ambiguous (0) or unambiguous (1)
alpha = 0.05; % Significance level
ambiguousIdxs = find(delaySig >= alpha);
fprintf('Indexes of ambiguous delays: %s\n', num2str(ambiguousIdxs));
save(boot_fname, 'delaySig', 'sigDists');

% Visualize non-zero and statistically ambiguous delays
plotGPparams_dlag(res.estParams, binWidth, rGroups, ...
    'plotAcross', true, ...
    'plotWithin', false, ...
    'units', 'ms', ...
    'sig', delaySig, ...
    'alpha', alpha);

%% ==========================================================
% 3b) Construct bootstrapped confidence intervals for latent delays
%     and timescales.
%  ================================================================

alpha = 0.05; % Construct (1-alpha) confidence intervals
bootParams = bootstrapGPparams(DLAG_struct, ...
    res.estParams, ...
    binWidth, ...
    numBootstrap, ...
    'alpha', alpha, ...
    'parallelize', parallelize, ...
    'numWorkers', numWorkers, ...
    'segLength', Inf, ...
    'tolLL', 1e-4, ...
    'maxIters', 10);
save(boot_fname, 'bootParams', '-append');
plotBootstrapGPparams_dlag(res.estParams, bootParams, binWidth, rGroups,...
    'overlayParams', false);

%% Save
cd(session_path)
% save(sprintf('%s_run%g_dim%g_dlag_workspace.mat', work, runIdx, i))
save(sprintf('%s_new_run%g_dim%g_bootstrap_%g.mat', work, runIdx, xDim_across, numBootstrap))
