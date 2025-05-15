function [bounds, reach_EMG_all, reach_all, up_thresh_all, peak_mag_all,...
    peak_all, ramp_all, trial_std_all] = ...
    get_reach_bounds(date,num_spouts,emg,success_times)
    % NEW
    % inputs the EMG (excluding Deltoid and Pectoralis) and success time
    % pairs (time or reward, time of touch), looks at 500ms before time of
    % reward and averages together the activity, when it increases a
    % certain amount of standard deviations from baseline, that is the 
    % start of the reach.
    
    %A prep time of 500ms is added to the start of the EMG trace from the
    %reward time in order to determine and EMG baseline and std.
    
    bounds = cell(1,num_spouts);
    ramp_all = cell(1,num_spouts);
    total_emg_avg = mean(emg,1);
    
    new_std = getGlobalSTD(date,total_emg_avg);
    %disp('Global Std. Dev. =')
    %disp(new_std)
    
    num_trials = 0;
    num_spouts_nonempty = [];
    for ii = 1:num_spouts
        if ~isempty(success_times{ii})
            num_trials = num_trials + length(success_times{ii}(:,1));
            num_spouts_nonempty = [num_spouts_nonempty,ii];
        end
    end
    
    prep = 500;
    reach_EMG_all = cell(num_trials,1);
    emg_post_rw = cell(num_trials,1);

    % Find the 30th percentile of EMG signal between reward and touch\
    trial_index = 0;
    for ii = num_spouts_nonempty
        for jj=1:size(success_times{ii},1)
            startt = success_times{ii}(jj,1);
            endt = success_times{ii}(jj,2);
            chunk = emg(:,startt-prep:endt);   
            % We are looking at the average muscles activity and so we take
            % the mean of the EMG signal across the 4 muscles.
            trace = mean(chunk,1);
            trial_index = trial_index+1;
            reach_EMG_all{trial_index} = trace;
        end
    end
    
    peak_all = zeros(1,num_trials);
    upper_thresh_array = [];
    trial_index = 0;
    for ii = 1:num_trials
        trial_index = trial_index+1;
        emg_post_rw{ii} = reach_EMG_all{ii}(501:end);
        upper_thresh_array = [upper_thresh_array,emg_post_rw{ii}]; 
        [~,peak_all(trial_index)] = max(emg_post_rw{ii});
        peak_all(trial_index) = peak_all(trial_index) + 500;
    end
    
    figure;
    histogram(upper_thresh_array);
    title('All EMG signals between reward and touch (all trials)');
    
    p_cutoff = 30;
    upper_thresh = prctile(upper_thresh_array,p_cutoff);
    %disp('Upper Threshold =')
    %disp(upper_thresh)
        
    reach_all = zeros(1,num_trials);
    up_thresh_all = zeros(1,num_trials);
    peak_mag_all = zeros(1,num_trials);
    trial_base_all = zeros(1,num_trials);
    trial_std_all = zeros(1,num_trials);
    trial_index = 0;
    for ii = num_spouts_nonempty  
        for jj=1:size(success_times{ii},1)
            trial_index = trial_index+1;
            startt = success_times{ii}(jj,1);
            endt = success_times{ii}(jj,2);
            full_chunk = emg(:,startt-prep:endt);   
            
            % We are looking at the average muscles activity and so we take
            % the mean of the EMG signal across the 4 muscles.
            full_trace = mean(full_chunk,1);
            baseline_trace = full_trace(1:prep-100);
            new_baseline = mean(baseline_trace);
            trial_base_all(trial_index) = new_baseline;
            
            trial_std = std(baseline_trace);
            trial_std_all(trial_index) = trial_std;

            [reach,up,ramp_time] = doubleThreshold(full_trace,...
                new_baseline,new_std,upper_thresh);
            if isempty(reach)
                reach = length(full_trace)-1;
            end
            
            [peak_mag_all(trial_index),~] = max(full_trace);
            up_thresh_all(trial_index) = up;
            reach_all(trial_index) = reach;
            ramp_all{ii} = [ramp_all{ii},ramp_time];
            
            % reach_start is the beginning of muscle movement and reach_end
            % is the time of the beam break.
            reach_start = startt-prep+reach;
            reach_end = endt;
            bound = [reach_start reach_end];
            bounds{ii} = [bounds{ii};bound];

        end    
    end
 end
    
    
