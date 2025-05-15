%%
clear all;
clc;
parpool('local', str2num(getenv('SLURM_NPROCS'))); %for QUEST ONLY

%cd '/home/dsb2139/ccm_for_quest'

addpath(genpath('.'))

%some parameters to set
run = 2
num_pairs = 10000;
num_nullmaps = 300;   

tau=1;
E=4;

tp_lags = -30:3:0;

num_trials_keep = 30;
pre = -200;
post = 799;

num_tps = length(tp_lags);

top10k_pairs = load('./10k_data/top_ten_thou_upd.mat').top_ten_thou_upd;
top10k_pairs = top10k_pairs(1:num_pairs, :);

c2r_coeffs = zeros(num_pairs, num_tps);
r2c_coeffs = zeros(num_pairs, num_tps);

best_c2r_coeffs = zeros(num_pairs, 1);
best_r2c_coeffs = zeros(num_pairs, 1);

best_c2r_tp = zeros(num_pairs, 1);
best_r2c_tp = zeros(num_pairs, 1);

c2r_null_coeffs = zeros(num_pairs, num_nullmaps);
r2c_null_coeffs = zeros(num_pairs, num_nullmaps);

smallest_num_trials = 70; 
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

        [CFA_FR, RFA_FR, reaches] = grab_pair(top10k_pairs(i,:));
    
        reaches = sort(reaches(random_trial_gen(i,:), :));%keep this many reaches
        total_samples = size(CFA_FR, 2);
        [bounds, reaching_bouts] = get_reaching_bouts(reaches, total_samples, pre, post);
    
        CFA_trial = fr_to_trial(CFA_FR, bounds);
        RFA_trial = fr_to_trial(RFA_FR, bounds);

        num_trials = length(CFA_trial);
        
        full_c2r_predic = [];
        full_c_orig = [];
    
        full_r2c_predic = [];
        full_r_orig = [];

        for j=1:num_trials
            [temp_c2r_predic, temp_c_orig] = ccm_qoneway(CFA_trial{1, j}, RFA_trial{1, j}, tau, E, tp_lags(k));
            [temp_r2c_predic, temp_r_orig] = ccm_qoneway(RFA_trial{1, j}, CFA_trial{1, j}, tau, E, tp_lags(k));

            full_c2r_predic = [full_c2r_predic; temp_c2r_predic];
            full_c_orig = [full_c_orig; temp_c_orig];

            full_r2c_predic = [full_r2c_predic; temp_r2c_predic];
            full_r_orig = [full_r_orig; temp_r_orig];
        end

        c2r_corr_temp=corrcoef(full_c_orig,full_c2r_predic,'Rows','complete');
        c2r_coeffs(i, k)=c2r_corr_temp(1,2);

        r2c_corr_temp=corrcoef(full_r_orig,full_r2c_predic,'Rows','complete');
        r2c_coeffs(i, k)=r2c_corr_temp(1,2);

    end
end
t=toc(t0);
fprintf('\n observed CCM took %f seconds', t)

save(sprintf('checkpoints/obs.mat'))

[best_c2r_coeffs, best_c2r_idx] = max(c2r_coeffs, [], 2);
[best_r2c_coeffs, best_r2c_idx] = max(r2c_coeffs, [], 2);

best_c2r_tp = tp_lags(best_c2r_idx);
best_r2c_tp = tp_lags(best_r2c_idx);

%fprintf
t0 = tic;
pb = progressbar('runinng CCM on null permutations...  ');
for k=1:num_nullmaps
    pb.print(k, num_nullmaps)
    parfor i=1:num_pairs
        [CFA_FR, RFA_FR, reaches] = grab_pair(top10k_pairs(i,:));
    
        reaches = sort(reaches(random_trial_gen(i,:), :));%keep this many reaches
        total_samples = size(CFA_FR, 2);
        
        [bounds, reaching_bouts] = get_reaching_bouts(reaches, total_samples, pre, post);

        [shift_CFA_trial, shift_RFA_trial] = circshift_pair(CFA_FR, RFA_FR, reaching_bouts);
        
        null_c2r_predic = [];
        null_c_orig = [];
    
        null_r2c_predic = [];
        null_r_orig = [];

        num_trials = length(shift_CFA_trial);
    
        for j=1:num_trials
            [nulltemp_c2r_predic, nulltemp_c_orig] = ccm_qoneway(shift_CFA_trial{1, j}, ...
                shift_RFA_trial{1, j}, tau, E, best_c2r_tp(i));
            [nulltemp_r2c_predic, nulltemp_r_orig] = ccm_qoneway(shift_RFA_trial{1, j}, ...
                shift_CFA_trial{1, j}, tau, E, best_r2c_tp(i));
    
            null_c2r_predic = [null_c2r_predic; nulltemp_c2r_predic];
            null_c_orig = [null_c_orig; nulltemp_c_orig];
            null_r2c_predic = [null_r2c_predic; nulltemp_r2c_predic];
            null_r_orig = [null_r_orig; nulltemp_r_orig];
        end
    
        c2r_nullcorr_temp=corrcoef(null_c_orig,null_c2r_predic,'Rows','complete');
        c2r_null_coeffs(i, k)=c2r_nullcorr_temp(1,2);
    
        r2c_nullcorr_temp=corrcoef(null_r_orig,null_r2c_predic,'Rows','complete');
        r2c_null_coeffs(i, k)=r2c_nullcorr_temp(1,2);

    end
    save(sprintf('checkpoints/nullcheckpoint_%d.mat', k))
end

t=toc(t0);
fprintf('null permutations took: %f seconsd', t)


[c2r_pvalues,r2c_pvalues] = get_CCM_pvalues(best_c2r_coeffs,best_r2c_coeffs, c2r_null_coeffs, r2c_null_coeffs);

save(sprintf('ccm__run%d_%dpairs_%dnullmaps.mat', run, num_pairs, num_nullmaps))

