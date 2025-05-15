clear all
close all

cd('Z:\Sarah');
addpath(genpath('DLAG'));

path = 'Z:\\Sarah\\DLAG';
run = 4;
dim = 4;
delay_max = 60; % 60 or 199

string = sprintf('CFA_RFA_struct_new_run%g_dim%g_bootstrap_1000', run, dim);


filelist = dir(fullfile(path, '**\*.*'));
session_list = [];


for k = 1:length(filelist)
    thisdir = filelist(k).name;
    if contains(thisdir, string) &&  contains(thisdir, '111')
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
ambig_RFA_var = [];
ambig_CFA_var = [];
nonzero_RFA_var = [];
nonzero_CFA_var = [];

CFA_across_var = {};
RFA_across_var = {};
CFA_within_var = {};
RFA_within_var = {};

CFA_across_vars_filt = {};
RFA_across_vars_filt = {};
CFA_within_vars_filt = {};
RFA_within_vars_filt = {};

figure;
hold on;
title(sprintf('Dimensionality %g', dim))
for i = 1:size(session,2)
    gp_params = plotGPparams_dlag_overlap(session{i}.res.estParams, session{i}.binWidth, session{i}.rGroups, ...
        'plotAcross', false, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', session{i}.delaySig, ...
        'alpha', session{i}.alpha);

    % Get delays
    delay = gp_params.DelayMatrix(session{i}.res.rGroups(2),:) ...
        - gp_params.DelayMatrix(session{i}.res.rGroups(1),:);

    ambig_idx = intersect(find(session{1,i}.delaySig >= 0.05), find(abs(delay) < 199));
    nonzero_idx = intersect(find(session{1,i}.delaySig < 0.05), find(abs(delay) < 199));
    
    temp = [delay(ambig_idx).', vars{i}{1}(ambig_idx).'];
    ambig_CFA_var = [ambig_CFA_var; temp];

    temp = [delay(nonzero_idx).', vars{i}{1}(nonzero_idx).'];
    nonzero_CFA_var = [nonzero_CFA_var; temp];

    temp = [delay(ambig_idx).', vars{i}{2}(ambig_idx).'];
    ambig_RFA_var = [ambig_RFA_var; temp];

    temp = [delay(nonzero_idx).', vars{i}{2}(nonzero_idx).'];
    nonzero_RFA_var = [nonzero_RFA_var; temp];

    % scatter_RFA_var(1+(i-1)*4:dim+(i-1)*4,1) = delay;
    % scatter_CFA_var(1+(i-1)*4:dim+(i-1)*4,1) = delay;
    % 
    % size(delay)

    % Indices of positive delays and negative delays
    % pos_idx = find(delay > 0);
    % neg_idx = find(delay < 0);
    % 
    % scatter_RFA_var(1+(i-1)*4:dim+(i-1)*4,2) = vars{i}{2}(1:dim);
    % scatter_CFA_var(1+(i-1)*4:dim+(i-1)*4,2) = vars{i}{1}(1:dim);
    % 
    CFA_across_var{i} = vars{i}{1}(1:dim);
    RFA_across_var{i} = vars{i}{2}(1:dim);

    CFA_within_var{i} = vars{i}{1}(end-run:end);
    RFA_within_var{i} = vars{i}{2}(end-run:end);


    CFA_across_var_filt{i} = CFA_across_var{i}((abs(delay) < delay_max)); 
    CFA_within_var_filt{i} = CFA_within_var{i}((abs(delay) < delay_max));
    RFA_across_var_filt{i} = RFA_across_var{i}((abs(delay) < delay_max));
    RFA_within_var_filt{i} = RFA_within_var{i}((abs(delay) < delay_max));

end

hold off;


%% Variance capture


figure;
hold on;
title(sprintf('Dimensionality %g', dim))
scatter(ambig_RFA_var(:,1), ambig_RFA_var(:,2), 'ro');
scatter(ambig_CFA_var(:,1), ambig_CFA_var(:,2), 'bo');
scatter(nonzero_RFA_var(:,1), nonzero_RFA_var(:,2), 'o', 'MarkerFaceColor','r');
scatter(nonzero_CFA_var(:,1), nonzero_CFA_var(:,2), 'o', 'MarkerFaceColor','b');
xline(0, '--');
xlim([-200, 200]);
xlabel('Delay (ms)');
ylabel('Proportion of shared variance captured');
legend('','','RFA variance', 'CFA variance');
hold off



%%
within_var_sums_CFA = [];
within_var_sums_RFA = [];
across_var_sums_CFA = [];
across_var_sums_RFA = [];

% All variance captured by modes
for i = 1:length(CFA_across_var_filt)
    within_var_sums_CFA(i,:) = sum(CFA_within_var_filt{1,i});
    within_var_sums_RFA(i,:) = sum(RFA_within_var_filt{1,i});
    across_var_sums_CFA(i,:) = sum(CFA_across_var_filt{1,i});
    across_var_sums_RFA(i,:) = sum(RFA_across_var_filt{1,i});
end


