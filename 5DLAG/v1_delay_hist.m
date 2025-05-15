clear all
close all

cd('Z:\Sarah');
addpath(genpath('DLAG'));

path = 'Z:\\Sarah\\DLAG';
run = 4;
dim = 4;
delay_max = 60; % 60 or 199

string = sprintf('CFA_RFA_struct_new_run%g_dim%g_workspace', run, dim);


filelist = dir(fullfile(path, '**\*.*'));
session_list = [];


for k = 1:length(filelist)
    thisdir = filelist(k).name;
    if contains(thisdir, string)
        session_list = [session_list; thisdir];
    end
end

session_list = unique(session_list, 'rows');

session = {};

for i = 1:size(session_list,1)
    session{i} = load(session_list(i,:));
end


%%

% Visualize non-zero and statistically ambiguous delays


cd('Z:\Sarah\DLAG');
addpath(genpath('DLAG-1.0.0'));

% Variance
vars = {};

for j = 1:length(session)
    vars{j} = session{j}.varexp.indiv;
end



% Plots
nonzeros = [];
ambig = [];
CC_var = {};
CR_var = {};
RC_var = {};
RR_var = {};
C2R_vars = {};
R2C_vars = {};
weights = {};
% CC_vars_filt = {};
% CR_vars_filt = {};
% RC_vars_filt = {};
% RR_vars_filt = {};
C2R_vars_filt = {};
R2C_vars_filt = {};

figure;
hold on;
title(sprintf('Dimensionality %g', dim))
for i = 1:size(session,2)
    gp_params = plotGPparams_dlag_overlap(session{i}.res.estParams, session{i}.binWidth, session{i}.rGroups, ...
        'plotAcross', true, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', session{i}.delaySig, ...
        'alpha', session{i}.alpha);

    % Get delays
    delay = gp_params.DelayMatrix(session{i}.res.rGroups(2),:) ...
        - gp_params.DelayMatrix(session{i}.res.rGroups(1),:);

    % Indices of positive delays and negative delays
    pos_idx = find(delay > 0);
    neg_idx = find(delay < 0);


    % Finding statistically significant delays
    ambig{i} = delay(session{1,i}.delaySig >= 0.05);
    nonzero_idx = find(session{1,i}.delaySig < 0.05);
    nonzeros{i} = delay(nonzero_idx);

    % Indices of sig. delays in specific directions
    pos_nonzero = intersect(pos_idx, nonzero_idx);
    neg_nonzero = intersect(neg_idx, nonzero_idx);

    % Postsynaptic variances of statistically significant delays
    C2R_vars{i} = vars{i}{2}(pos_nonzero);
    R2C_vars{i} = vars{i}{1}(neg_nonzero);

    weights{i}(pos_nonzero) = C2R_vars{i};
    weights{i}(neg_nonzero) = R2C_vars{i};

    C2R_vars_filt{i} = C2R_vars{i}((abs(delay(pos_nonzero)) < delay_max)); 
    R2C_vars_filt{i} = R2C_vars{i}((abs(delay(neg_nonzero)) < delay_max));


end


%% Removing 200 ms

all_nonzeros = [];
all_weights = [];


for i = 1:size(nonzeros,2)
    temp = nonzeros{i};
    temp2 = weights{i};
    all_nonzeros = [all_nonzeros temp];
    all_weights = [all_weights temp2];
    all_weights = all_weights(find(all_weights));
    nonzeros_filt = all_nonzeros(abs(all_nonzeros) < delay_max);
    weights_filt = all_weights(abs(all_nonzeros) < delay_max);
end

med_nonzero = median(nonzeros_filt);

% edges = -200:20:200;
edges = -100:20:100;
[histw, ~] = histwc(nonzeros_filt, ones(1, length(nonzeros_filt)), edges);


% Delays
figure;
hold on
title(sprintf('Dimensionality %g', dim))
b = bar(edges, histw, 'histc');
b.FaceColor = [0.9,0.6,0.6];
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
xlim([-100, 100])
hold off


% Weighted

[histw, ~] = histwc(nonzeros_filt, weights_filt, edges);

figure;
hold on
title(sprintf('Dimensionality %g', dim))
b = bar(edges, histw, 'histc');
b.FaceColor = [0.9,0.6,0.6];
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
ylabel('Weighted counts')
xlim([-100, 100])
hold off


%% Variance capture

C2R_var_sums = [];
R2C_var_sums = [];

% All variance captured by modes with significant delays
for i = 1:length(C2R_vars_filt)
    C2R_var_sums(i,:) = sum(C2R_vars_filt{1,i});
    R2C_var_sums(i,:) = sum(R2C_vars_filt{1,i});
end


C2R_mean_var = mean(C2R_var_sums);
R2C_mean_var = mean(R2C_var_sums);

C2R_sem = std(C2R_var_sums)/sqrt(numel(C2R_var_sums));
R2C_sem = std(R2C_var_sums)/sqrt(numel(R2C_var_sums));

% Plot
figure;
hold on
errorbar(1, R2C_mean_var, R2C_sem, 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r', 'Color', 'r');
errorbar(2, C2R_mean_var, C2R_sem, 'o', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'b', 'Color', 'b');
xlim([0, 3]);
ylim([0, 1]);
ylabel('Proportion of shared variance captured');
legend('RFA to CFA', 'CFA to RFA');
hold off


save_path = uigetdir('', 'Save Location');
cd(save_path)
save(sprintf('run%g_dim%g_long_new_var_stats.mat', run, dim), 'run', 'dim','C2R_mean_var', 'R2C_mean_var', 'C2R_sem', 'R2C_sem');

