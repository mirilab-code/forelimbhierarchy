clear all
close all


path = 'Z:\\Sarah\\multiprobe_DLAG';
dim = 4;
num_sessions = 3; %will change

session{1} = load(sprintf('%s\\mp27_03282023\\mp27_03282023_CFA_RFA_pcs_run3_no_cv_dim%g_workspace.mat', path, dim));
session{2} = load(sprintf('%s\\mp29_04042023\\mp29_04042023_CFA_RFA_pcs_run3_no_cv_dim%g_workspace.mat', path, dim));
session{3} = load(sprintf('%s\\mp29_04052023\\mp29_04052023_CFA_RFA_pcs_run3_no_cv_dim%g_workspace.mat', path, dim));

save('delay_hist_session_load.mat')

%%

% Visualize non-zero and statistically ambiguous delays

% Plots
nonzeros = [];
ambig = [];

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
    
    delay = gp_params.DelayMatrix(session{i}.res.rGroups(2),:) ...
           - gp_params.DelayMatrix(session{i}.res.rGroups(1),:);
    ambig{i} = delay(session{1,i}.delaySig >= 0.05);
    nonzeros{i} = delay(session{1,i}.delaySig < 0.05);
end


%% Removing 200 ms 

all_nonzeros = [];
for i = 1:size(nonzeros,2)
    temp = nonzeros{i};
    all_nonzeros = [all_nonzeros temp];
    nonzeros_filt = all_nonzeros(abs(all_nonzeros) < 199);
end
med_nonzero = median(nonzeros_filt);

figure;
hold on
title(sprintf('Dimensionality %g', dim))
histogram(nonzeros_filt, 20)
xline(med_nonzero, '--r')
