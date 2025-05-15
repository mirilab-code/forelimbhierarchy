
clear all;

%% Load files
session_path = uigetdir('', 'Data Path');
cd(session_path)
date = session_path(end-7:end);

spike_trains = uigetfile('','Spike Trains');
load(spike_trains);

reach_bounds = uigetfile('','Reach Bounds');
load(reach_bounds);


%% Binary spike matrix

pre = 200;
post = 200;

CFA_trains = pyr_CFA_RFA_trains_sort{1};
RFA_trains = pyr_CFA_RFA_trains_sort{2};

CFA_binary = logical([]);
RFA_binary = logical([]);


% Find endpoint
CFA_end_times = [];
for i = 1:size(CFA_trains,1)
    CFA_end_times(i,:) = round(CFA_trains{i}(end));
end
last1 = max(CFA_end_times);

RFA_end_times = [];
for i = 1:size(RFA_trains,1)
    RFA_end_times(i,:) = round(RFA_trains{i}(end));
end
last2 = max(RFA_end_times);

last = max(last1,last2);

% CFA
for n = 1:size(CFA_trains,1)
    trn = CFA_trains{n};
    binary = logical(zeros(1, last));
    if round(trn(1)) == 0
        trn(1) = 1;
    end
    binary(round(trn)) = 1;
    CFA_binary = [CFA_binary; binary];
end



% RFA
for n = 1:size(RFA_trains,1)
    trn = RFA_trains{n};
    binary = logical(zeros(1, last));
    if round(trn(1)) == 0
        trn(1) = 1;
    end
    binary(round(trn)) = 1;
    RFA_binary = [RFA_binary; binary];
end

% full_binary = [CFA_binary; RFA_binary];

%% FR exclusions

CFA_bad = zeros(size(CFA_binary,1),1);
RFA_bad = zeros(size(RFA_binary,1),1);


% CFA exclusions
for x = 1:size(CFA_binary, 1)
    current_neuron = CFA_binary(x,:);
    CFA_trials = logical([]);
    for j = 1:size(reach_bounds_edit, 1)
        CFA_windows = current_neuron(:, reach_bounds_edit(j)-pre:reach_bounds_edit(j)+post);
        CFA_trials = [CFA_trials; CFA_windows];
    end
    BadTrials = zeros(size(CFA_trials, 1), 1);
    BadTrials(find(sum(CFA_trials,2) == 0)) = 1;
    if sum(BadTrials) > 0
        CFA_bad(x) = 1;
    end
end


% RFA exclusions
for x = 1:size(RFA_binary, 1)
    current_neuron = RFA_binary(x,:);
    RFA_trials = logical([]);
    for j = 1:size(reach_bounds_edit, 1)
        RFA_windows = current_neuron(:, reach_bounds_edit(j)-pre:reach_bounds_edit(j)+post);
        RFA_trials = [RFA_trials; RFA_windows];
    end
    BadTrials = zeros(size(RFA_trials, 1), 1);
    BadTrials(find(sum(RFA_trials,2) == 0)) = 1;
    if sum(BadTrials) > 0
        RFA_bad(x) = 1;
    end
end


CFA_excluded = CFA_binary(find(CFA_bad == 0),:);
RFA_excluded = RFA_binary(find(RFA_bad == 0),:);


% CFA_exc = logical([]);
% RFA_exc = logical([]);
% 
% for j = 1:size(reach_bounds_edit,1)
%     onset_ms = reach_bounds_edit(j);
%     CFA_windows = CFA_binary(:,onset_ms-pre:onset_ms+post);
%     RFA_windows = RFA_binary(:,onset_ms-pre:onset_ms+post);
%     CFA_exc = [CFA_exc CFA_windows];
%     RFA_exc = [RFA_exc RFA_windows];
% end
% 
% figure
% histogram(sum(CFA_exc,2)/length(CFA_exc), 100)
% figure
% histogram(sum(RFA_exc,2)/length(RFA_exc), 100)
% 
% CFA_excluded = CFA_binary(sum(CFA_exc,2)/length(CFA_exc) > 0.001,:);
% RFA_excluded = RFA_binary(sum(RFA_exc,2)/length(RFA_exc) > 0.001,:);
% 
% 
full_binary = [CFA_excluded; RFA_excluded];


%% Struct creation


all_struct = struct;
for t = 1:size(reach_bounds_edit,1)
    all_struct(t).trialId = t;
    onset_ms = reach_bounds_edit(t);
    onset_windows = full_binary(:,onset_ms-pre:onset_ms+post);
    all_struct(t).spikes = onset_windows;
end



%% getSeq

cd('Z:\Sarah\DLAG\')
addpath('DLAG-1.0.0\')

all_seq = getSeq(all_struct, 1);


%% Save
save_path = uigetdir('', 'Save Location');
cd(save_path)

save(sprintf('%s_binned_spike_counts.mat', date), 'all_seq', 'CFA_excluded', 'RFA_excluded');


%% Notes

% Need: xDims_across and xDims_within from FA?

