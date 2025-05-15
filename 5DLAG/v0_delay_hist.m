clear all
close all

cd('Z:\Sarah');
addpath(genpath('DLAG'));

path = 'Z:\\Sarah\\DLAG';
run = 4;
dim = 4;

string = sprintf('CFA_RFA_struct_run%g_dim%g_workspace', run, dim);


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
nonzeros_scaled = [];
ambig = [];
C2R_vars = {};
R2C_vars = {};
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

    % Scaling delays by postsynaptic variance
%  fix   delay_scaled(pos_idx) = delay(pos_idx).*vars{i}{2}(pos_idx); 
%     delay_scaled(neg_idx) = delay(neg_idx).*vars{i}{1}(neg_idx);

    % Finding statistically significant delays
    ambig{i} = delay(session{1,i}.delaySig >= 0.05);
    nonzero_idx = find(session{1,i}.delaySig < 0.05);
    nonzeros{i} = delay(nonzero_idx);
    nonzeros_scaled{i} = delay_scaled(nonzero_idx);

    % Indices of sig. delays in specific directions
    pos_nonzero = intersect(pos_idx, nonzero_idx);
    neg_nonzero = intersect(neg_idx, nonzero_idx);

    % Postsynaptic variances of statistically significant delays
    C2R_vars{i} = vars{i}{2}(pos_nonzero);
    R2C_vars{i} = vars{i}{1}(neg_nonzero);


    C2R_vars_filt{i} = C2R_vars{i}((abs(delay(pos_nonzero)) < 199));
    R2C_vars_filt{i} = R2C_vars{i}((abs(delay(neg_nonzero)) < 199));

end


%% Removing 200 ms

all_nonzeros = [];
all_nonzeros_scaled = [];


for i = 1:size(nonzeros,2)
    temp = nonzeros{i};
    temp2 = nonzeros_scaled{i};
    all_nonzeros = [all_nonzeros temp];
    all_nonzeros_scaled = [all_nonzeros_scaled temp2];
    nonzeros_filt = all_nonzeros(abs(all_nonzeros) < 199);
    nonzeros_scaled_filt = all_nonzeros_scaled(abs(all_nonzeros) < 199);
end
med_nonzero = median(nonzeros_filt);

edges = -200:20:200;

% Delays
figure;
hold on
title(sprintf('Dimensionality %g', dim))
histogram(nonzeros_filt, edges, 'FaceColor', [0.9,0.6,0.6])
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
xlim([-200, 200])
hold off

% Scaled by shared var
edges2 = -100:10:100;

figure;
hold on
title(sprintf('Dimensionality %g', dim))
histogram(nonzeros_scaled_filt, edges2, 'FaceColor', [0.9,0.6,0.6])
xline(median(nonzeros_scaled_filt), '--r')
xlabel('Delay scaled (ms)')
xlim([-100, 100])
hold off

%% Variance capture
%{
CFA_vars = [];
RFA_vars = [];

C2R_vars{i} = vars{i}{2}(pos_nonzero);
R2C_vars{i} = vars{i}{1}(neg_nonzero);

% Variances captured by each area
for i = 1:length(vars)
    CFA_vars(i,:) = vars{i}{1, 1};
    RFA_vars(i,:) = vars{i}{1, 2};
end
%}
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

length(nonzeros_filt)/(18*dim)

save_path = uigetdir('', 'Save Location');
cd(save_path)
save(sprintf('run%g_dim%g_var_stats.mat', run, dim), 'run', 'dim','C2R_mean_var', 'R2C_mean_var', 'C2R_sem', 'R2C_sem');

