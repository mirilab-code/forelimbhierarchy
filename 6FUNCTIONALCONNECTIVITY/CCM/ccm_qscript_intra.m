%%
clear all;
clc;
%parpool('local', str2num(getenv('SLURM_NPROCS'))); %for QUEST ONLY


addpath(genpath('.'))
pairpath = './10k_data/intra_pairs';

%some parameters to set
run = 0;
num_pairs = 100;
num_nullmaps = 50;   

tau=1;
E=4;

tp_lags = -30:3:0;

num_trials_keep = 30;
pre = -200;
post = 500;

num_tps = length(tp_lags);

intra_pairs = load(pairpath).intra_pairs;
intra_pairs = intra_pairs(1:num_pairs, :);

x2y_coeffs = zeros(num_pairs, num_tps);
y2x_coeffs = zeros(num_pairs, num_tps);

x2y_null_coeffs = zeros(num_pairs, num_nullmaps);
y2x_null_coeffs = zeros(num_pairs, num_nullmaps);

smallest_num_trials = 70; %this is quite stupid, but lets not go down a rabbit hole.
random_trial_gen = zeros(num_pairs, num_trials_keep);
for i=1:num_pairs
    random_trial_temp = randperm(smallest_num_trials);
    random_trial_gen(i, :) = random_trial_temp(1:num_trials_keep);
end

t0=tic;

pb = progressbar('running CCM on different time lags...  ');
for k=1:num_tps
    pb.print(k, num_tps);
    tp = tp_lags(k);
    parfor i=1:num_pairs

        [x_FR, y_FR, reaches] = grab_pair_intra(intra_pairs(i,:));
    
        reaches = sort(reaches(random_trial_gen(i,:), :));%keep this many reaches
        total_samples = size(x_FR, 2);
        [bounds, reaching_bouts] = get_reaching_bouts(reaches, total_samples, pre, post);
    
        x_trial = fr_to_trial(x_FR, bounds);
        y_trial = fr_to_trial(y_FR, bounds);

        num_trials = length(x_trial);
        
        full_x2y_predic = [];
        full_x_orig = [];
    
        full_y2x_predic = [];
        full_y_orig = [];

        for j=1:num_trials
            [temp_x2y_predic, temp_x_orig] = ccm_qoneway(x_trial{1, j}, y_trial{1, j}, tau, E, tp_lags(k));
            [temp_y2x_predic, temp_y_orig] = ccm_qoneway(y_trial{1, j}, x_trial{1, j}, tau, E, tp_lags(k));

            full_x2y_predic = [full_x2y_predic; temp_x2y_predic];
            full_x_orig = [full_x_orig; temp_x_orig];

            full_y2x_predic = [full_y2x_predic; temp_y2x_predic];
            full_y_orig = [full_y_orig; temp_y_orig];
        end

        x2y_corr_temp=corrcoef(full_x_orig,full_x2y_predic,'Rows','complete');
        x2y_coeffs(i, k)=x2y_corr_temp(1,2);

        y2x_corr_temp=corrcoef(full_y_orig,full_y2x_predic,'Rows','complete');
        y2x_coeffs(i, k)=y2x_corr_temp(1,2);

    end
end
t=toc(t0);
fprintf('\n observed CCM took %f seconds', t)

%save(sprintf('checkpoints/obs.mat'))

[best_x2y_coeffs, best_x2y_idx] = max(x2y_coeffs, [], 2);
[best_y2x_coeffs, best_y2x_idx] = max(y2x_coeffs, [], 2);

best_x2y_tp = tp_lags(best_x2y_idx);
best_y2x_tp = tp_lags(best_y2x_idx);

%fprintf
t0 = tic;
pb = progressbar('runinng CCM on null permutations...  ');
for k=1:num_nullmaps
    pb.print(k, num_nullmaps)
    parfor i=1:num_pairs
        [x_FR, y_FR, reaches] = grab_pair_intra(intra_pairs(i,:));
    
        reaches = sort(reaches(random_trial_gen(i,:), :));%keep this many reaches
        total_samples = size(x_FR, 2);
        
        [bounds, reaching_bouts] = get_reaching_bouts(reaches, total_samples, pre, post);

        [shift_x_trial, shift_y_trial] = circshift_pair(x_FR, y_FR, reaching_bouts);
        
        null_x2y_predic = [];
        null_x_orig = [];
    
        null_y2x_predic = [];
        null_y_orig = [];

        num_trials = length(shift_x_trial);
    
        for j=1:num_trials
            [nulltemp_x2y_predic, nulltemp_x_orig] = ccm_qoneway(shift_x_trial{1, j}, ...
                shift_y_trial{1, j}, tau, E, best_x2y_tp(i));
            [nulltemp_y2x_predic, nulltemp_y_orig] = ccm_qoneway(shift_y_trial{1, j}, ...
                shift_x_trial{1, j}, tau, E, best_y2x_tp(i));
    
            null_x2y_predic = [null_x2y_predic; nulltemp_x2y_predic];
            null_x_orig = [null_x_orig; nulltemp_x_orig];
            null_y2x_predic = [null_y2x_predic; nulltemp_y2x_predic];
            null_y_orig = [null_y_orig; nulltemp_y_orig];
        end
    
        x2y_nullcorr_temp=corrcoef(null_x_orig,null_x2y_predic,'Rows','complete');
        x2y_null_coeffs(i, k)=x2y_nullcorr_temp(1,2);
    
        y2x_nullcorr_temp=corrcoef(null_y_orig,null_y2x_predic,'Rows','complete');
        y2x_null_coeffs(i, k)=y2x_nullcorr_temp(1,2);

    end
    %save(sprintf('checkpoints/nullcheckpoint_%d.mat', k))
end

t=toc(t0);
fprintf('null permutations took: %f seconsd', t)


[y2x_pvalues,x2y_pvalues] = get_CCM_pvalues(best_x2y_coeffs,best_y2x_coeffs, x2y_null_coeffs, y2x_null_coeffs, best_x2y_tp, best_y2x_tp);

%save(sprintf('ccm__run%d_%dpairs_%dnullmaps.mat', run, num_pairs, num_nullmaps))

