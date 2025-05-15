clear all
close all

session_path = uigetdir('', 'Data Path');
cd(session_path)
date = session_path(end-7:end);


dims = [2,4,6,8];

for i = dims
    dim{dims == i} = load(sprintf('%s_CFA_RFA_pcs_run3_no_cv_dim%g_workspace.mat', date, i));
end


%%

% Visualize non-zero and statistically ambiguous delays

% Separate plots
% for i = 1:8
%     plotGPparams_dlag(dim{i}.res.estParams, dim{i}.binWidth, dim{i}.rGroups, ...
%         'plotAcross', true, ...
%         'plotWithin', false, ...
%         'units', 'ms', ...
%         'sig', dim{i}.delaySig, ...
%         'alpha', dim{i}.alpha);
% end


% Together

% Colors
% mycolors = {[138 43 226]./255,... % Purple
%     [30 144 255]./255,... % Blue
%     [64 224 208]./255,... % Aqua
%     [0 255 0]./255,... % Green
%     [255 215 0]./255,... % Yellow
%     [255 140 0]./255,... % Orange
%     [255 0 0]./255,... % Red
%     [255 20 147]./255}; %Pink


mycolors = {[138 43 226]./255,... % Purple, dim2
    [64 224 208]./255,... % Aqua, dim4
    [255 215 0]./255,... % Yellow, dim6
    [255 0 0]./255,... % Red, dim8
    [255 20 147]./255}; %Pink

% Plots
nonzeros = [];
ambig = [];

figure;
hold on;
title(sprintf('Session %s', date))
for i = 1:size(dim,2)
    gp_params = plotGPparams_dlag_overlap(dim{i}.res.estParams, dim{i}.binWidth, dim{i}.rGroups, ...
        'plotAcross', true, ...
        'plotWithin', false, ...
        'units', 'ms', ...
        'sig', dim{i}.delaySig, ...
        'alpha', dim{i}.alpha, ...
        'color', mycolors{i});
    
    delay = gp_params.DelayMatrix(2,:);
    ambig{i} = delay(dim{1,i}.delaySig >= 0.05);
    nonzeros{i} = delay(dim{1,i}.delaySig < 0.05);
end


%{
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
%}
%% Removing 200 ms 

all_nonzeros = [];
for i = 1:size(nonzeros,2)
    temp = nonzeros{i};
    all_nonzeros = [all_nonzeros temp];
    nonzeros_filt = all_nonzeros(abs(all_nonzeros) < 199);
end
