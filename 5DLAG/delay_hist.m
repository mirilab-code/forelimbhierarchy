clear all
close all

cd('Z:\Sarah');
addpath(genpath('DLAG'));

path = 'Z:\\Sarah\\DLAG';
run = 4;
dim = 4;
delay_max = 60; % 60 or 199
animal_nums = [1,2,3,3,4,4,4,5,5,5,5,6,6,6,6];

string = sprintf('CFA_RFA_struct_new_run%g_dim%g_bootstrap', run, dim);


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
median_delays = [];

CC_all = {};
CR_all = {};
RC_all = {};
RR_all = {};

CC_var = {}; % CFA lead, CFA var
CR_var = {}; % CFA lead, RFA var
RC_var = {}; % RFA lead, CFA var
RR_var = {}; % RFA lead, RFA var

CFA_within_var = {};
RFA_within_var = {};

CFA_weights = {};
RFA_weights = {};

CC_vars_filt = {};
CR_vars_filt = {};
RC_vars_filt = {};
RR_vars_filt = {};



figure;
hold on;
title(sprintf('Dimensionality %g', dim))
for i = 1:size(session,2)

    animal = animal_nums(i);

    CFA_within_var{animal,i} = session{1,i}.varexp.indiv{1,1}(end-dim+1:end);
    RFA_within_var{animal,i} = session{1,i}.varexp.indiv{1,2}(end-dim+1:end);

    % CFA_within_var{i} = session{1,i}.varexp.indiv{1,1}(end-dim+1:end);
    % RFA_within_var{i} = session{1,i}.varexp.indiv{1,2}(end-dim+1:end);

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

    median_delays(i) = median(nonzeros{i}(abs(nonzeros{i}) < delay_max));

    % Indices of sig. delays in specific directions
    pos_nonzero = intersect(pos_idx, nonzero_idx);
    neg_nonzero = intersect(neg_idx, nonzero_idx);

    % Postsynaptic variances of all delays
    CR_all{i} = vars{i}{2}(pos_nonzero);
    RC_all{i} = vars{i}{1}(neg_nonzero);

    % Presynaptic variances of all delays
    CC_all{i} = vars{i}{1}(pos_nonzero);
    RR_all{i} = vars{i}{2}(neg_nonzero);


    % Postsynaptic variances of statistically significant delays
    CR_var{i} = vars{i}{2}(pos_nonzero);
    RC_var{i} = vars{i}{1}(neg_nonzero);

    % Presynaptic variances of statistically significant delays
    CC_var{i} = vars{i}{1}(pos_nonzero);
    RR_var{i} = vars{i}{2}(neg_nonzero);

    CFA_weights{i}(pos_nonzero) = CC_var{i}; % CFA leading, CFA variance 
    CFA_weights{i}(neg_nonzero) = RC_var{i}; % RFA leading, CFA variance

    RFA_weights{i}(pos_nonzero) = CR_var{i}; % CFA leading, RFA variance
    RFA_weights{i}(neg_nonzero) = RR_var{i}; % RFA leading, RFA variance


    CR_vars_filt{animal,i} = CR_var{i}((abs(delay(pos_nonzero)) < delay_max)); 
    RC_vars_filt{animal,i} = RC_var{i}((abs(delay(neg_nonzero)) < delay_max));
    CC_vars_filt{animal,i} = CC_var{i}((abs(delay(pos_nonzero)) < delay_max)); 
    RR_vars_filt{animal,i} = RR_var{i}((abs(delay(neg_nonzero)) < delay_max));

    % CR_vars_filt{i} = CR_var{i}((abs(delay(pos_nonzero)) < delay_max)); 
    % RC_vars_filt{i} = RC_var{i}((abs(delay(neg_nonzero)) < delay_max));
    % CC_vars_filt{i} = CC_var{i}((abs(delay(pos_nonzero)) < delay_max)); 
    % RR_vars_filt{i} = RR_var{i}((abs(delay(neg_nonzero)) < delay_max));


end


%% Removing 200 ms

all_nonzeros = [];
all_weights_CFA = [];
all_weights_RFA = [];


for i = 1:size(nonzeros,2)
    temp = nonzeros{i};
    temp2 = CFA_weights{i};
    temp3 = RFA_weights{i};
    all_nonzeros = [all_nonzeros temp];
    all_weights_CFA = [all_weights_CFA temp2];
    all_weights_RFA = [all_weights_RFA temp3];

    all_weights_CFA = all_weights_CFA(find(all_weights_CFA));
    all_weights_RFA = all_weights_RFA(find(all_weights_RFA));
    nonzeros_filt = all_nonzeros(abs(all_nonzeros) < delay_max);
    weights_filt_CFA = all_weights_CFA(abs(all_nonzeros) < delay_max);
    weights_filt_RFA = all_weights_RFA(abs(all_nonzeros) < delay_max);
end

med_nonzero = median(nonzeros_filt);
[p, h, stats] = signrank(median_delays, 0, 'tail', 'left');

mean_CFA_weight = sum(nonzeros_filt.*weights_filt_CFA)/sum(weights_filt_CFA);
mean_RFA_weight = sum(nonzeros_filt.*weights_filt_RFA)/sum(weights_filt_RFA);


edges = -200:20:200;
% edges = -100:20:100;
[histw, ~] = histwc(nonzeros_filt, ones(1, length(nonzeros_filt)), edges);


% Delays
figure;
hold on
title(sprintf('Dimensionality %g', dim))
b = bar(edges, histw, 'histc');
b.FaceColor = [0.9,0.6,0.6];
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
xlim([-200, 200])
hold off


% Weighted by CFA var

[histw, ~] = histwc(nonzeros_filt, weights_filt_CFA, edges);

figure;
hold on
title(sprintf('Dimensionality %g', dim))
b = bar(edges, histw, 'histc');
b.FaceColor = [0.9,0.6,0.6];
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
ylabel('Weighted counts (CFA)')
% xlim([-100, 100])
xlim([-200, 200])
hold off

% Weighted by RFA var

[histw, ~] = histwc(nonzeros_filt, weights_filt_RFA, edges);

figure;
hold on
title(sprintf('Dimensionality %g', dim))
b = bar(edges, histw, 'histc');
b.FaceColor = [0.9,0.6,0.6];
xline(med_nonzero, '--r')
xlabel('Delay (ms)')
ylabel('Weighted counts (RFA)')
% xlim([-100, 100])
xlim([-200, 200])
hold off


%% Variance capture

CR_var_sums = [];
RC_var_sums = [];
RR_var_sums = [];
CC_var_sums = [];
RFA_var_sums = [];
CFA_var_sums = [];

% All variance captured by modes with significant delays

for i = 1:size(CR_vars_filt, 1)
CR = [];
RC = [];
CC = [];
RR = [];
CFA = [];
RFA = [];

    for j = 1:size(CR_vars_filt, 2)
        CR(:,j) = sum(CR_vars_filt{i,j});
        RC(:,j) = sum(RC_vars_filt{i,j});
        CC(:,j) = sum(CC_vars_filt{i,j});
        RR(:,j) = sum(RR_vars_filt{i,j});
        CFA(:,j) = sum(CFA_within_var{i,j});
        RFA(:,j) = sum(RFA_within_var{i,j});
    end
    CR_var_sums(i,:) = sum(CR)/nnz(CR);
    RC_var_sums(i,:) = sum(RC)/nnz(RC);
    CC_var_sums(i,:) = sum(CC)/nnz(CC);
    RR_var_sums(i,:) = sum(RR)/nnz(RR);
    CFA_var_sums(i,:) = sum(CFA)/nnz(CFA);
    RFA_var_sums(i,:) = sum(RFA)/nnz(RFA);
end


CR_var_sums(isnan(CR_var_sums)) = 0;
RC_var_sums(isnan(RC_var_sums)) = 0;
CC_var_sums(isnan(CC_var_sums)) = 0;
RR_var_sums(isnan(RR_var_sums)) = 0;
CFA_var_sums(isnan(CFA_var_sums)) = 0;
RFA_var_sums(isnan(RFA_var_sums)) = 0;

% for i = 1:length(CR_vars_filt)
%     CR_var_sums(i,:) = sum(CR_vars_filt{1,i});
%     RC_var_sums(i,:) = sum(RC_vars_filt{1,i});
%     CC_var_sums(i,:) = sum(CC_vars_filt{1,i});
%     RR_var_sums(i,:) = sum(RR_vars_filt{1,i});
% end

% for i = 1:size(CFA_within_var,2)
%     CFA_var_sums(i,:) = sum(CFA_within_var{1,i});
%     RFA_var_sums(i,:) = sum(RFA_within_var{1,i});
% end

RFA_data = [CR_var_sums, RR_var_sums, RFA_var_sums];
CFA_data = [CC_var_sums, RC_var_sums, CFA_var_sums];

% RFA_data = {CR_var_sums, RR_var_sums, RFA_var_sums};
% CFA_data = {CC_var_sums, RC_var_sums, CFA_var_sums};

% RFA var plot
figure;
hold on
plot([1,2,3], RFA_data, 'o-k', 'MarkerFaceColor','k');
xlim([0, 4]);
xticks([1,2,3]);
xticklabels({'CFA->RFA', 'RFA->CFA', 'Within RFA'});
ylim([-0.05, 1]);
ylabel('Summed prop. of shared variance capture');
title('RFA')
hold off

% figure;
% hold on
% boxplotGroup(RFA_data, 'primaryLabels', {'CFA->RFA', 'RFA->CFA', 'Within RFA'}, 'symbol', '');
% xlim([0, 4]);
% ylim([-0.05, 1]);
% ylabel('Summed prop. of shared variance capture');
% title('RFA');
% hold off


% CFA var plot
figure;
hold on
plot([1,2,3], CFA_data, 'o-k', 'MarkerFaceColor','k');
xlim([0, 4]);
xticks([1,2,3]);
xticklabels({'CFA->RFA', 'RFA->CFA', 'Within CFA'});
ylim([-0.05, 1]);
ylabel('Summed prop. of shared variance capture');
title('CFA')
hold off


% figure;
% hold on
% boxplotGroup(CFA_data, 'primaryLabels', {'CFA->RFA', 'RFA->CFA', 'Within CFA'}, 'symbol', '');
% xlim([0, 4]);
% ylim([-0.05, 1]);
% ylabel('Summed prop. of shared variance capture');
% title('CFA')
% hold off


