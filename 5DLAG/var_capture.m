clear all
close all

cd('Z:\Sarah');
addpath(genpath('DLAG'));

path = 'Z:\\Sarah\\DLAG\\overall_new';

% Within = dim optimal, across = run optimal
x = 'within'; % across or within to go on x-axis, the other will be set to optimal

string = 'long_new_var_stats.mat';
optimal = 4;

filelist = dir(fullfile(path, '**\*.*'));
results_list = [];


for k = 1:length(filelist)
    thisdir = filelist(k).name;
    if strcmp(x, 'across') == 1
        if contains(thisdir, string) && contains(thisdir, sprintf('run%g', optimal))
            results_list = [results_list; thisdir];
        end
    else
        if contains(thisdir, string) && contains(thisdir, sprintf('dim%g', optimal))
            results_list = [results_list; thisdir];
        end
    end
end

results_list = unique(results_list, 'rows');


results = {};

for i = 1:size(results_list,1)
    results{i} = load(results_list(i,:));
end

%% Plot
C2R_stats = [];
R2C_stats = [];

for i = 1:size(results, 2)
    if strcmp(x, 'across') == 1
        C2R_stats(i,1) = results{i}.dim;
        R2C_stats(i,1) = results{i}.dim;
    else
        C2R_stats(i,1) = results{i}.run;
        R2C_stats(i,1) = results{i}.run;
    end
    C2R_stats(i,2) = results{i}.C2R_mean_var;
    C2R_stats(i,3) = results{i}.C2R_sem;

    R2C_stats(i,2) = results{i}.R2C_mean_var;
    R2C_stats(i,3) = results{i}.R2C_sem;
end

figure;
hold on
errorbar(R2C_stats(:,1), R2C_stats(:,2), R2C_stats(:,3), 'o', 'MarkerFaceColor','r', 'MarkerEdgeColor', 'r', 'Color', 'r');
errorbar(C2R_stats(:,1), C2R_stats(:,2), C2R_stats(:,3), 'o', 'MarkerFaceColor','b', 'MarkerEdgeColor', 'b', 'Color', 'b');
xlim([1.5, 6.5]);
ylim([0, 0.5]);
ylabel('Proportion of shared variance captured');
if strcmp(x, 'across')
    xlabel('Across-area dimensionalities')
    title(sprintf('Within-area dimensionality %g', optimal))
else
    xlabel('Within-area dimensionalities')
    title(sprintf('Across-area dimensionality %g', optimal))
end
legend('RFA to CFA', 'CFA to RFA');
hold off


