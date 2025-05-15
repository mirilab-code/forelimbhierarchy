clear all
close all




for i = 1:4
    session_path = uigetdir('', sprintf('Data Path %g', i));
    cd(session_path)
    date = session_path(end-7:end);

    dim{i} = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim4_workspace.mat', date));
    dates(i,:) = date;
end



% dim1 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim1_workspace.mat', date));
% dim2 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim2_workspace.mat', date));
% dim3 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim3_workspace.mat', date));
% dim4 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim4_workspace.mat', date));
% dim5 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim5_workspace.mat', date));
% dim6 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim6_workspace.mat', date));
% dim7 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim7_workspace.mat', date));
% dim8 = load(sprintf('%s_CFA_first_pcs_run3_no_cv_dim8_workspace.mat', date));

%%

% Visualize non-zero and statistically ambiguous delays

% Separate plots
for i = 1:size(dim, 2)
    plotGPparams_dlag_fixed_lim(dim{i}.res.estParams, dim{i}.binWidth, dim{i}.rGroups, ...
        'plotAcross', true, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', dim{i}.delaySig, ...
        'alpha', dim{i}.alpha, ...
        'session', dates(i,:));
end


%% Together

mycolors = {[138 43 226]./255,... % Purple
    [30 144 255]./255,... % Blue
    [0 255 0]./255,... % Green
    [255 140 0]./255,... % Orange
    [255 0 0]./255}; % Red


figure;
hold on;
title(sprintf('Session %s', date))
for i = 1:size(dim,2)
    plotGPparams_dlag_overlap(dim{i}.res.estParams, dim{i}.binWidth, dim{i}.rGroups, ...
        'plotAcross', true, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', dim{i}.delaySig, ...
        'alpha', dim{i}.alpha, ...
        'color', mycolors{i});
end



% Ambiguous delays
figure;
hold on;
title(sprintf('Session %s', date))
for i = 1:size(dim,2)
    plotGPparams_dlag_ambig(dim{i}.res.estParams, dim{i}.binWidth, dim{i}.rGroups, ...
        'plotAcross', true, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', dim{i}.delaySig, ...
        'alpha', dim{i}.alpha, ...
        'color', mycolors{i});
end

