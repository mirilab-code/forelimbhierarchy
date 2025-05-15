function [all_trials,all_trials_averages] = FRtrialAverage(FR,trial_bounds)
%INPUT:
%   FR: a cell array containing the FR time series
%   trial_bounds: a num_trials X 2 matrix containing the start and end times
%   for each trial
%OUTPUT:
%   all_trials: a (num_units x trial_duration x total_num_trials) matrix 
%       containing the FR time series for each unit for every trial.
%   all_trials_averages: a (num_units x trial_duration) matrix containing
%       the  FR time series for each units averaged across trials.

% Trim the FR matrix to only the windows around the reach initiation
% Arranged in a (num_units x trial_duration x total_num_trials) matrix 
% containing the FR time series for each unit for every trial.
num_probes = length(FR);
all_trials = cell(num_probes,1);
for ii = 1:num_probes
    for jj = 1:length(FR{ii}(:,1))
        for kk = 1:length(trial_bounds(:,1))
            all_trials{ii}(jj,:,kk) = FR{ii}(jj,...
                trial_bounds(kk,1):trial_bounds(kk,2));
        end
    end
end

% Find the average FR for each unit across all trials.
all_trials_averages = cell(num_probes,1);
for ii = 1:num_probes
    all_trials_averages{ii} = mean(all_trials{ii},3);
end
end

