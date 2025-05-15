function [sigCFA, sigRFA] = find_neurmuscle_correlations(data_folder)

% disp('loading things...');
disp(data_folder);
files = fullfile(data_folder,'*.mat');
matFiles = dir(files);


% do neurons
events_ind = find(strcmp({matFiles.name},'events.mat'));
events_fname = fullfile(data_folder, matFiles(events_ind).name);
data = load(events_fname);
events = data.events;

% then emg
emg_ind = find(strcmp({matFiles.name},'EMG.mat'));
emg_fname = fullfile(data_folder, matFiles(emg_ind).name);
data = load(emg_fname);

EMG = data.EMG;
train_CFA = events_to_train(events{1});
train_CFA = train_CFA(~cellfun('isempty',train_CFA));
train_RFA = events_to_train(events{2});
train_RFA = train_RFA(~cellfun('isempty',train_RFA));

disp('done loading')
%%
duration = size(EMG,2);
FRcfa = trains_to_firingrate(train_CFA,duration);
FRrfa = trains_to_firingrate(train_RFA,duration);

%%
each_muscle_corr_CFA = [];
each_muscle_corr_RFA = [];
for i=1:size(EMG,1)
    this_muscle = EMG(i,:);
    stackCFA = [this_muscle; FRcfa];
    corr_cfa = corr(stackCFA');

    stackRFA = [this_muscle; FRrfa];
    corr_rfa = corr(stackRFA');
    
    this_neur_corr_CFA = corr_cfa(1,2:end);
    this_neur_corr_RFA = corr_rfa(1,2:end);

    each_muscle_corr_CFA = cat(1,each_muscle_corr_CFA,this_neur_corr_CFA);
    each_muscle_corr_RFA = cat(1,each_muscle_corr_RFA,this_neur_corr_RFA);
end


%% now do the bootstrapping
% do just one half for some reason
% ht = round(duration/2);
% EMG = EMG(:,1:ht);
% FRcfa = FRcfa(:,1:ht);
% FRrfa = FRrfa(:,1:ht);
%

niter = 300;
bootstrap_CFA = bootstrap_corr(EMG,FRcfa,niter,data_folder);
bootstrap_RFA = bootstrap_corr(EMG,FRrfa,niter,data_folder);


% make sure to do a two-tailed test to find correlated and anticorrelated neurons
%% doing steps in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3624014/#:~:text=Calculate%20the%20serial,and%20anti%2D%20correlation.

% do for each muscle
muscles_pvals_CFA = [];
muscles_pvals_RFA = [];
for i=1:size(EMG,1)
    this_bootstrap_CFA = squeeze(bootstrap_CFA(1,:,:));
    this_bootstrap_RFA = squeeze(bootstrap_RFA(1,:,:));
    this_obs_CFA = each_muscle_corr_CFA(i,:);
    this_obs_RFA = each_muscle_corr_RFA(i,:);

    % fraction of bootstraps larger than observed values
    frac_CFA = mean(this_obs_CFA' < this_bootstrap_CFA,2);
    frac_RFA = mean(this_obs_RFA' < this_bootstrap_RFA,2);

    % turn those into p values
    pCFA = 2* (0.5 - abs(mean(this_bootstrap_CFA < frac_CFA,2) -0.5))';
    pRFA = 2* (0.5 - abs(mean(this_bootstrap_RFA < frac_RFA,2) -0.5))';

    muscles_pvals_CFA = [muscles_pvals_CFA; pCFA];
    muscles_pvals_RFA = [muscles_pvals_RFA; pRFA];

end

% neurons that are significantly correlated with *at least one* muscle are
% marked as significant
sigCFA = any(muscles_pvals_CFA<0.05,1);
sigRFA = any(muscles_pvals_RFA<0.05,1);

%%
sname = split(data_folder,'\');
sname = ['muscleneurocorrelations\' sname{end-1} '.mat'];

save(sname,'EMG','sigCFA','sigRFA','train_CFA','train_RFA','bootstrap_CFA','bootstrap_RFA','niter');
fprintf('saved to %s \n',sname);




end