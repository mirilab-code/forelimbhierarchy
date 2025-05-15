addpath(genpath('Z:\Scripts'));

% Input 'RFA', 'CFA', or 'ALL'
animal_subset = 'ALL'

if strcmp(animal_subset,'RFA')
    animals = {'re7','re8','re9'};
elseif strcmp(animal_subset,'CFA')
    animals = {'re10','re12','re14'};
elseif strcmp(animal_subset,'ALL')
    animals = {'re7','re8','re9','re10','re12','re14'};
elseif strcmp(animal_subset,'NEW')
    animals = {'re7','re9'};
else
    animals = {animal_subset};
end

num_animals = length(animals);
abs_diff_emg_agg = [];
musc = 4;

for animal_index = 1:num_animals

    if strcmp(animals{animal_index},'re7')
        dates = ['10132021';'10142021';'10152021';'10182021';'10192021'];
    elseif strcmp(animals{animal_index},'re8')
        dates = ['10272021';'10282021';'10292021';'11012021'];
    elseif strcmp(animals{animal_index},'re9')
        dates = ['11162021';'11172021';'11182021';'11192021'];
    elseif strcmp(animals{animal_index},'re10')
        dates = ['01182022';'01192022';'01202022';'01212022'];
    elseif strcmp(animals{animal_index},'re12')
        dates = ['12162021';'12172021'];
    elseif strcmp(animals{animal_index},'re14')
        dates = ['03012021';'03022022';'03032022';'03042022';'03072022'];
    else
        disp('Error! Wrong date!')
    end

    num_dates = size(dates,1);
    EMG_all_control = [];
    EMG_all_stim = [];
    fr_all_control = [];
    fr_all_stim = [];
    control_lag = [];
    control_pre = [];
    stim_lag = [];
    stim_pre = [];
    for date_index = 1:num_dates
        animal = animals{animal_index};
        date = dates(date_index,:);
        str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\EMG.mat'];
        load(str)
        str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\analogin.mat'];
        load(str)
        str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\digitalin.mat'];
        load(str)
        if strcmp('111192021',date)
            str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\neurons_curated.mat'];
        else
            str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\neurons.mat'];
        end
        load(str)
        str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\events.mat'];
        load(str)
        
        num_muscles = size(EMG,1);
        for ii = 1:num_muscles
            EMG(ii,:) = zscore(EMG(ii,:));
        end
        
        if ~strcmp('111192021',date)
            N = neurons{1};
        end
        d  = [N.depth];
        exclude = d>2000;
        N(exclude) = [];
        w = [N.width];
        exclude = w<14;
        %N(exclude) = [];
        trains = {N.train};
        depths = [N.depth];
        w = [N.width];

        fr = neurons_to_firingrate_oneSided(N);
%         train = events_to_train(events{1});
%         fr = train_to_firingrate(train);
%         good_units = zeros(length(train),1);
%         for ii = 1:length(train)
%             good_units(ii) = ~isempty(train{ii});
%         end
%         fr = fr(logical(good_units),1:size(EMG,2));
        num_neurons = size(fr,1);
        full_duration = size(fr,2);

        bin_spike = zeros(num_neurons,full_duration);
        for ii = 1:num_neurons
            temp = round(trains{ii});
            temp(temp==0) = [];
            bin_spike(ii,temp) = 1;
        end

        touch_sensor1 = 4;
        touch_sensor2 = 3;
        touch_sensor3 = 2;
        touch_sensor4 = 1;
        resting_bar_sensor = 5;
        reward1 = 6;
        reward2 = 7;
        reward3 = 8;
        reward4 = 9;
        suction = 10;
        % Create reward cell containing the times when the reward sound is
        % detected (this is for both successful and unsuccesful trials).
        rw = get_reward_times(digitalin, [reward1 reward2 reward3 reward4]);
    
        % correct for the reward sound by subtracting 100
        rw1 = rw{reward1}-100;
        rw2 = rw{reward2}-100;
        rw3 = rw{reward3}-100;
        rw4 = rw{reward4}-100;
        rw_corr = {rw1,rw2,rw3,rw4};
        num_spouts = 4;
        touch_channel1 = digitalin(touch_sensor1,:);
        touch_channel2 = digitalin(touch_sensor2,:);
        touch_channel3 = digitalin(touch_sensor3,:);
        touch_channel4 = digitalin(touch_sensor4,:);
        touch_channels = [touch_channel1; touch_channel2; ...
            touch_channel3; touch_channel4];
        touch_chan_add = sum(touch_channels,1);
        touch_chan_add(touch_chan_add > 0) = 1;
        allowable_error = 50;
        success_bounds = find_success(num_spouts,rw_corr,touch_channels,...
            digitalin(suction,:),allowable_error);
        all_success_nospout = [];
        for ii = 1:num_spouts
            all_success_nospout = [all_success_nospout; success_bounds{ii}];
        end

        summed_EMG = sum(EMG,1);
        
        stim_channel = 2;
        control_channel = 7;
        stim_trace = analogin(stim_channel,:);
        control_trace = analogin(control_channel,:);
        stim_trace = stim_trace>1.5;
        control_trace = control_trace>1.5;
        stim_times = find(diff(stim_trace)>0);
        control_times = find(diff(control_trace)>0);
        good_control_inds = [];
        good_stim_inds = [];
        
        for ii = 1:size(all_success_nospout,1)
            tone = all_success_nospout(ii,1);
            touch = all_success_nospout(ii,2);
            for jj = 1:length(control_times)
                check = control_times(jj);
                if check > tone && check < touch
                    good_control_inds = [good_control_inds,jj];
                    control_lag = [control_lag,check-tone];
                    control_pre = [control_pre,touch-check];
                end
            end
            for jj = 1:length(stim_times)
                check = stim_times(jj);
                if check > tone && check < touch
                    good_stim_inds = [good_stim_inds,jj];
                    stim_lag = [stim_lag,check-tone];
                    stim_pre = [stim_pre,touch-check];
                end
            end
        end
        

%         stim_times = stim_times(good_stim_inds);
%         control_times = control_times(good_control_inds);
            
        window = [-500,500];
        time = window(1):window(2);
        duration = length(time);
        
        bad_inds = [];
        for ii = 1:length(stim_times)
            if stim_times(ii)+window(1)<0 || stim_times(ii)+window(2)>size(EMG,2)
                bad_inds = [bad_inds,ii];
            end
        end
        if ~isempty(bad_inds)
            stim_times(bad_inds) = [];
        end
        
        bad_inds = [];
        for ii = 1:length(control_times)
            if control_times(ii)+window(1)<0 || control_times(ii)+window(2)>size(EMG,2)
                bad_inds = [bad_inds,ii];
            end
        end
        if ~isempty(bad_inds)
            control_times(bad_inds) = [];
        end
               
        % figure;
        % plot(stim_trace)
        % hold on
        % plot(control_trace)
        
        EMG_stim = zeros(length(stim_times),duration,num_muscles);
        EMG_control = zeros(length(control_times),duration,num_muscles);
        for ii = 1:length(stim_times)
            for jj = 1:num_muscles
                EMG_stim(ii,:,jj) = ...
                    EMG(jj,stim_times(ii)+window(1):stim_times(ii)+window(2));
            end
        end
        for ii = 1:length(control_times)
            for jj = 1:num_muscles
                EMG_control(ii,:,jj) = ...
                    EMG(jj,control_times(ii)+window(1):control_times(ii)+window(2));
            end
        end
        EMG_all_control = cat(1,EMG_all_control,EMG_control);
        EMG_all_stim = cat(1,EMG_all_stim,EMG_stim);
        
        bin_stim = zeros(length(stim_times),duration,num_neurons);
        bin_control = zeros(length(control_times),duration,num_neurons);
        fr_stim = zeros(length(stim_times),duration,num_neurons);
        fr_control = zeros(length(control_times),duration,num_neurons);
        for ii = 1:length(stim_times)
            for jj = 1:num_neurons
                fr_stim(ii,:,jj) = ...
                    fr(jj,stim_times(ii)+window(1):stim_times(ii)+window(2));
                bin_stim(ii,:,jj) = ...
                    bin_spike(jj,stim_times(ii)+window(1):stim_times(ii)+window(2));
            end
        end
        for ii = 1:length(control_times)
            for jj = 1:num_neurons
                fr_control(ii,:,jj) = ...
                    fr(jj,control_times(ii)+window(1):control_times(ii)+window(2));
                bin_control(ii,:,jj) = ...
                    bin_spike(jj,control_times(ii)+window(1):control_times(ii)+window(2));
            end
        end

        % for ii = 1:size(EMG_control,1)
        %     figure;
        %     for jj = 1:num_muscles
        %         plot(EMG_control(ii,:,jj))
        %         hold on
        %     end
        %     plot(sum(EMG_control(ii,:,:),3),'k--')
        % end
        % 
        % for ii = 1:size(EMG_stim,1)
        %     figure;
        %     for jj = 1:num_muscles
        %         plot(EMG_stim(ii,:,jj))
        %         hold on
        %     end
        %     plot(sum(EMG_stim(ii,:,:),3),'c--')
        % end
        
        % figure;
        % for ii = 1:size(EMG_control,1)
        %     plot(sum(EMG_control(ii,:,:),3),'k')
        %     hold on
        % end
        % figure;
        % for ii = 1:size(EMG_stim,1)
        %     plot(sum(EMG_stim(ii,:,:),3),'c')
        %     hold on
        % end
%         
%         figure;
%         for ii = 1:size(fr_control,1)
%             figure
%             neuron_mean = mean(fr_control(ii,:,:),3);
%             plot(-100:100,neuron_mean(400:600),'k')
%             hold on
%         end
%         figure;
%         for ii = 1:size(fr_stim,1)
%             figure
%             neuron_mean = mean(fr_stim(ii,:,:),3);
%             plot(-100:100,neuron_mean(400:600),'c')
%             hold on
%         end


        EMG_control_sum = sum(EMG_control,3);
        EMG_stim_sum = sum(EMG_stim,3);
        
        SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
        SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));
        
        EMG_control_trial_avg = mean(EMG_control_sum);
        EMG_stim_trial_avg = mean(EMG_stim_sum);
        EMG_control_bl_mean = mean(EMG_control_trial_avg(1:300));
        EMG_stim_bl_mean = mean(EMG_stim_trial_avg(1:300));
        EMG_control_sub = EMG_control_trial_avg - EMG_control_bl_mean;
        EMG_stim_sub = EMG_stim_trial_avg - EMG_stim_bl_mean;

        fr_control_mean = mean(fr_control,3);
        fr_stim_mean = mean(fr_stim,3);

        fr_all_control = cat(1,fr_all_control,fr_control_mean);
        fr_all_stim = cat(1,fr_all_stim,fr_stim_mean);
        
        fr_SEM_control = std(fr_control_mean,1)/sqrt(size(fr_control_mean,1));
        fr_SEM_stim = std(fr_stim_mean,1)/sqrt(size(fr_stim_mean,1));
        
        fr_control_trial_avg = mean(fr_control_mean);
        fr_stim_trial_avg = mean(fr_stim_mean);
        fr_control_bl_mean = mean(fr_control_trial_avg(1:300));
        fr_stim_bl_mean = mean(fr_stim_trial_avg(1:300));
        fr_control_sub = fr_control_trial_avg - fr_control_bl_mean;
        fr_stim_sub = fr_stim_trial_avg - fr_stim_bl_mean;

        separations = zeros(1,duration);
        for ii = 1:duration
            separations(ii) = (EMG_control_sub(:,ii)-SEM_control(ii))-...
                (EMG_stim_sub(:,ii)+SEM_stim(ii));
        end
        first_diverge = find(separations(502:end)>0,1,'first');
        disp(first_diverge)
        
        %swin = 1:1001;
        swin = 401:601;
        figure;
        boundedline(time(swin),EMG_control_sub(swin),SEM_control(swin),'k','alpha')
        hold on
        boundedline(time(swin),EMG_stim_sub(swin),SEM_stim(swin),'c','alpha')
        str = [animal,' - ',date,': Trial Averaged EMG Activty Surrounding Laser Inactivation'];
        title(str)
        xlabel('Inactivation at t = 0 ms')
        ylabel('EMG Actvity')
        str1 = ['n = ',int2str(length(control_times)), ' trials'];
        str2 = ['n = ',int2str(length(stim_times)), ' trials'];
        legend(str1,'control',str2,'stim','Location','southeast')

        figure;
        boundedline(time(swin),fr_control_sub(swin),fr_SEM_control(swin),'k','alpha')
        hold on
        boundedline(time(swin),fr_stim_sub(swin),fr_SEM_stim(swin),'c','alpha')
        str = [animal,' - ',date,': Trial Averaged CFA Activty Surrounding Laser Inactivation'];
        title(str)
        xlabel('Inactivation at t = 0 ms')
        ylabel('Neural Actvity')
        str1 = ['n = ',int2str(length(control_times)), ' trials'];
        str2 = ['n = ',int2str(length(stim_times)), ' trials'];
        legend(str1,'control',str2,'stim','Location','southeast')

   
    end


%% Individual and Aggregate Muscle Activity
EMG_control_agg_mouse = [];
EMG_stim_agg_mouse = [];

EMG_mouse_control = {};
EMG_mouse_stim = {};

for musc = 1:num_muscles
    EMG_control_sum = EMG_all_control(:,:,musc);
    EMG_stim_sum = EMG_all_stim(:,:,musc);
    
    SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
    SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));
    
    EMG_control_trial_avg = mean(EMG_control_sum);
    EMG_stim_trial_avg = mean(EMG_stim_sum);
    EMG_control_bl_mean = 0 %mean(EMG_control_trial_avg(1:300));
    EMG_stim_bl_mean = 0 %mean(EMG_stim_trial_avg(1:300));
    EMG_control_sub = EMG_control_trial_avg - EMG_control_bl_mean;
    EMG_stim_sub = EMG_stim_trial_avg - EMG_stim_bl_mean;
    
    separations = zeros(1,duration);
    for ii = 1:duration
        separations(ii) = (EMG_control_sub(:,ii)-SEM_control(ii))-...
            (EMG_stim_sub(:,ii)+SEM_stim(ii));
    end
    first_diverge = find(separations(502:end)>0,1,'first');
    disp(first_diverge)
    
    % Aggregate muscle activity from all muscles
    EMG_control_agg_mouse(musc, :) = EMG_control_sub(swin);
    EMG_stim_agg_mouse(musc, :) = EMG_stim_sub(swin);
    SEM_control_agg_mouse(musc, :) = SEM_control(swin);
    SEM_stim_agg_mouse(musc, :) = SEM_stim(swin);

    EMG_mouse_control{1,musc} = EMG_control_sum;
    EMG_mouse_stim{1,musc} = EMG_stim_sum;

    %figure;
    subplot(1,4,musc);
    boundedline(time(swin),EMG_control_sub(swin),SEM_control(swin),'k','alpha')
    hold on
    boundedline(time(swin),EMG_stim_sub(swin),SEM_stim(swin),'c','alpha')
    muscles_names = {'Bicep','Tricep','ECR','PL'};
    str = [animal,' - ',muscles_names{musc}];%,': Trial Averaged EMG Activty Surrounding Inactivation'];
    title(str)
    ylim([-0.5 3]);
    xlabel('Inactivation at t = 0 ms')
    ylabel('EMG Actvity')
    str1 = ['n = ',int2str(size(EMG_control_sum,1)), ' trials'];
    str2 = ['n = ',int2str(size(EMG_stim_sum,1)), ' trials'];
%     legend(str1,'control',str2,'stim','Location','southeast')
end

% Aggregate muscle activity from all animals and muscles
EMG_control_agg{animal_index,:} = EMG_control_agg_mouse;
EMG_stim_agg{animal_index,:} = EMG_stim_agg_mouse;
SEM_control_agg{animal_index,:} = SEM_control_agg_mouse;
SEM_stim_agg{animal_index,:} = SEM_stim_agg_mouse;

EMG_control_all_trials(animal_index,:) = EMG_mouse_control;
EMG_stim_all_trials(animal_index,:) = EMG_mouse_stim;

%
EMG_control_sum = sum(EMG_all_control,3);
EMG_stim_sum = sum(EMG_all_stim,3);

SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));

EMG_control_trial_avg = mean(EMG_control_sum);
EMG_stim_trial_avg = mean(EMG_stim_sum);
EMG_control_bl_mean = 0 %mean(EMG_control_trial_avg(1:300));
EMG_stim_bl_mean = 0 %mean(EMG_stim_trial_avg(1:300));
EMG_control_sub = EMG_control_trial_avg - EMG_control_bl_mean;
EMG_stim_sub = EMG_stim_trial_avg - EMG_stim_bl_mean;

fr_control_mean = fr_all_control;
fr_stim_mean = fr_all_stim;

fr_SEM_control = std(fr_control_mean,1)/sqrt(size(fr_control_mean,1));
fr_SEM_stim = std(fr_stim_mean,1)/sqrt(size(fr_stim_mean,1));

fr_control_trial_avg = mean(fr_control_mean);
fr_stim_trial_avg = mean(fr_stim_mean);
fr_control_bl_mean = mean(fr_control_trial_avg(1:300));
fr_stim_bl_mean = mean(fr_stim_trial_avg(1:300));
fr_control_sub = fr_control_trial_avg - fr_control_bl_mean;
fr_stim_sub = fr_stim_trial_avg - fr_stim_bl_mean;

separations = zeros(1,duration);
for ii = 1:duration
    separations(ii) = (EMG_control_sub(:,ii)-SEM_control(ii))-...
        (EMG_stim_sub(:,ii)+SEM_stim(ii));
end
first_diverge = find(separations(502:end)>0,1,'first');
disp(first_diverge)

figure;
boundedline(time(swin),EMG_control_sub(swin),SEM_control(swin),'k','alpha')
hold on
boundedline(time(swin),EMG_stim_sub(swin),SEM_stim(swin),'c','alpha')
str = [animal,' - Combined Sessions: Trial Averaged EMG Activty Surrounding Laser Inactivation'];
title(str)
xlabel('Inactivation at t = 0 ms')
ylabel('EMG Actvity')
str1 = ['n = ',int2str(size(EMG_control_sum,1)), ' trials'];
str2 = ['n = ',int2str(size(EMG_stim_sum,1)), ' trials'];
legend(str1,'control',str2,'stim','Location','southeast')

EMG_control_agg_musc_avg{animal_index,:} = EMG_control_sub(swin);
EMG_stim_agg_musc_avg{animal_index,:} = EMG_stim_sub(swin);
SEM_control_agg_musc_avg{animal_index,:} = SEM_control(swin);
SEM_stim_agg_musc_avg{animal_index,:} = SEM_stim(swin);

figure;
boundedline(time(swin),fr_control_sub(swin),fr_SEM_control(swin),'k','alpha')
hold on
boundedline(time(swin),fr_stim_sub(swin),fr_SEM_stim(swin),'c','alpha')
str = [animal,' - Combined Sessions: Trial Averaged CFA Activty Surrounding Laser Inactivation'];
title(str)
xlabel('Inactivation at t = 0 ms')
ylabel('Neural Actvity')
str1 = ['n = ',int2str(size(fr_all_control,1)), ' trials'];
str2 = ['n = ',int2str(size(fr_all_stim,1)), ' trials'];
legend(str1,'control',str2,'stim','Location','southeast')

%%
% max_EMG = max(summed_EMG);
% figure;
% hold on
% plot(summed_EMG)
% plot(max_EMG*touch_chan_add)
% plot(max_EMG*stim_trace)
% plot(max_EMG*control_trace)
max_lag = max([stim_lag,control_lag]);
max_pre = max([stim_pre,control_pre]);
figure;
histogram(stim_lag,0:max_lag/40:max_lag)
hold on
histogram(control_lag,0:max_lag/40:max_lag)
legend('stim','control')
title('Laser-Tone')
figure;
histogram(stim_pre,0:max_pre/40:max_pre)
hold on
histogram(control_pre,0:max_pre/40:max_pre)
legend('stim','control')
title('Touch-Laser')
%% Absolute Difference
abs_diff_emg = abs(EMG_control_sub-EMG_stim_sub);
abs_diff_emg = abs_diff_emg - median(abs_diff_emg(401:450));
abs_diff_emg = abs_diff_emg ./ 4;
figure;
plot(time(swin),abs_diff_emg(swin))
title('Absolute Difference in EMG Activity (Control and Stim)')

abs_diff_cfa = abs(fr_control_sub-fr_stim_sub);
abs_diff_cfa = abs_diff_cfa - median(abs_diff_cfa(401:450));
figure;
hold on
plot(time(swin),abs_diff_cfa(swin))
title('Absolute Difference in CFA Activity (Control and Stim)')

abs_diff_emg_agg(animal_index,:) = abs_diff_emg(swin);

%%
% swin = 251:751;
% fr_control_pca = fr_control_mean(:,swin);
% fr_stim_pca = fr_stim_mean(:,swin);
% time = -(length(swin)-1)/2:(length(swin)-1)/2;
% 
% init_bl =100;
% [princomp_control,score_control,latent_control] = pca(fr_control_pca');
% score_control_mod = score_control(:,1) - mean(score_control(1:init_bl,1));
% score_control_mod = score_control_mod/max(score_control_mod);
% figure;
% plot(time,score_control_mod,'k')
% 
% [princomp_stim,score_stim,latent_stim] = pca(fr_stim_pca');
% score_stim_mod = score_stim(:,1) - mean(score_stim(1:init_bl,1));
% score_stim_mod = score_stim_mod/max(score_stim_mod);
% hold on
% plot(time,score_stim_mod,'c')
% xline(0)
% %xlim([time(1),time(end)])
% xlim([-100,100])
% xlabel('Time (ms)')
% str = [animal,': First CFA neural mode (trial averaged) surrounding inactivation'];
% title(str)
% legend('control','stim','Location','Southeast')

end

%% Load EMG Trials
clear all; close all;

cd('Z:\Sajishnu\MATLAB\Variables');
load('EMG_control_all_trials.mat');
load('EMG_stim_all_trials.mat');

addpath(genpath('Z:\Scripts'));

% Input 'RFA', 'CFA', or 'ALL'
animal_subset = 'ALL'

if strcmp(animal_subset,'RFA')
    animals = {'re7','re8','re9'};
elseif strcmp(animal_subset,'CFA')
    animals = {'re10','re12','re14'};
elseif strcmp(animal_subset,'ALL')
    animals = {'re7','re8','re9','re10','re12','re14'};
else
    animals = {animal_subset};
end

musc = 4;
num_animals = length(animals);
%% Trial Exclusion

EMG_control_all_trials_exclude = EMG_control_all_trials;
EMG_stim_all_trials_exclude = EMG_stim_all_trials;
% Reorganize Data to trials x -50-0ms per muscle
% for ii = 1:musc
%     for i = 1:num_animals
%         EMG_control_all_trials_exclude{i,ii} = reshape(EMG_control_all_trials_exclude{i,ii}, 1, size(EMG_control_all_trials_exclude{i,ii},2), size(EMG_control_all_trials_exclude{i,ii},1));
%         EMG_stim_all_trials_exclude{i,ii} = reshape(EMG_stim_all_trials_exclude{i,ii}, 1, size(EMG_stim_all_trials_exclude{i,ii},2), size(EMG_stim_all_trials_exclude{i,ii},1));
%     end
% end
% 
for i = 1:num_animals
    EMG_control_all_trials_exclude{i,1} = [EMG_control_all_trials_exclude{i,1}(:,451:501),EMG_control_all_trials_exclude{i,2}(:,451:501),EMG_control_all_trials_exclude{i,3}(:,451:501),EMG_control_all_trials_exclude{i,4}(:,451:501)];
    EMG_stim_all_trials_exclude{i,1} = [EMG_stim_all_trials_exclude{i,1}(:,451:501),EMG_stim_all_trials_exclude{i,2}(:,451:501),EMG_stim_all_trials_exclude{i,3}(:,451:501),EMG_stim_all_trials_exclude{i,4}(:,451:501)];
end

for i = 1:3
    EMG_control_all_trials_exclude(:,2) = [];
    EMG_stim_all_trials_exclude(:,2) = [];
end

% for i = 1:num_animals
%     EMG_control_all_trials_exclude{i,1} = mean(EMG_control_all_trials_exclude{i,1},1);
%     s = size(EMG_control_all_trials_exclude{i,1});
%     EMG_control_all_trials_exclude{i,1} = reshape(EMG_control_all_trials_exclude{i,1}, s(3), s(2));
% 
%     EMG_stim_all_trials_exclude{i,1} = mean(EMG_stim_all_trials_exclude{i,1},1);
%     s = size(EMG_stim_all_trials_exclude{i,1});
%     EMG_stim_all_trials_exclude{i,1} = reshape(EMG_stim_all_trials_exclude{i,1}, s(3), s(2));   
% end


% Concatenate control and stim
for i = 1:num_animals
    EMG_control_stim_all_trials_exclude{i,1} = cat(1, EMG_control_all_trials_exclude{i,1},EMG_stim_all_trials_exclude{i,1});
end

% Find outliers
thld = 1; % Multiplier of std
bad_trials_list = [];
distance_metric = 'euclidean';
for i = 1:num_animals
    % pdist2 funciton
    test = EMG_control_stim_all_trials_exclude{i,1};
    dist = pdist2(test, test, distance_metric);

    % mean resulting matrix
    result = mean(dist);
        
%         figure(i);
%         histogram(result)
%         hold on
%         title('mean distance of trials from other trials')

    % remove outlier trials
    threshold = mean(result)+(thld*std(result)); % what do I want to set this as?
    bad_trials = find(result>threshold);

    bad_trials_control = bad_trials(bad_trials < size(EMG_control_all_trials_exclude{i,1},1));
    bad_trials_stim = bad_trials(bad_trials > size(EMG_control_all_trials_exclude{i,1},1)) - size(EMG_control_all_trials_exclude{i,1},1);

    for ii = 1:musc
        EMG_control_all_trials{i,ii}(bad_trials_control,:) = [];
        EMG_stim_all_trials{i,ii}(bad_trials_stim,:) = [];
    end
    bad_trials_list(i,1) = length(bad_trials_stim);
    disp(['Animal ',num2str(i),' - ',num2str(length(bad_trials_stim)/size(EMG_stim_all_trials_exclude{i,1},1))])
end

disp('bad trials gone!')

%% Individual Muscles
muscles_names = {'Bicep','Tricep','ECR','PL'};
animals = {'re7','re8','re9','re10','re12','re14'};
swin = [401:601];
subplot_pos = 1;
figure();
sgtitle([distance_metric, num2str(thld),' * STD Surrounding the Mean '])
for ii = 1:musc
    for i = 1:num_animals
        % Get muscle trials
        EMG_control_sum = EMG_control_all_trials{i,ii};
        EMG_stim_sum = EMG_stim_all_trials{i,ii};
        % Calc SEM
        SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
        SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));
        % Avg trials
        EMG_control_avg = mean(EMG_control_sum,1);
        EMG_stim_avg = mean(EMG_stim_sum,1);
        % BL subtract
        EMG_control_avg = EMG_control_avg - mean(EMG_control_avg(401:476));
        EMG_stim_avg = EMG_stim_avg - mean(EMG_stim_avg(401:476));

        subplot(musc,num_animals,subplot_pos);
        boundedline(swin,EMG_control_avg(swin),SEM_control(swin),'k','alpha')
        hold on
        boundedline(swin,EMG_stim_avg(swin),SEM_stim(swin),'c','alpha')
        muscles_names = {'Bicep','Tricep','ECR','PL'};
        animals = {'re7','re8','re9','re10','re12','re14'};
        str = [animals{i},' - ',muscles_names{ii}, num2str(bad_trials_list(i,1)),'/',num2str(size(EMG_control_all_trials{i,1},1))];%,': Trial Averaged EMG Activty Surrounding Inactivation'];
        title(str)
        ylim([-0.5 3]);
        xlabel('Inactivation at t = 0 ms')
        ylabel('EMG Actvity')
        str1 = ['n = ',int2str(size(EMG_control_sum,1)), ' trials'];
        str2 = ['n = ',int2str(size(EMG_stim_sum,1)), ' trials'];
%         legend(str1,'control',str2,'stim','Location','southeast')
        hold on

        subplot_pos = subplot_pos + 1;
    end
end

%% Plot Individual Animals Muscles
% RFA Inac Used RE9
% CFA Inac Used RE14
musc_animal = 're14';
muscles_names = {'Bicep','Tricep','ECR','PL'};
animals = {'re7','re8','re9','re10','re12','re14'};
swin = [451:901];
subplot_pos = 1;
if strcmp(musc_animal,'re7')
    i = 1;
elseif strcmp(musc_animal,'re8')
    i = 2;
elseif strcmp(musc_animal,'re9')
    i = 3;
elseif strcmp(musc_animal,'re10')
    i = 4;
elseif strcmp(musc_animal,'re12')
    i = 5;
elseif strcmp(musc_animal,'re14')
    i = 6;
else
    disp('Error');
end
figure();
% sgtitle([animals{i},' - ',distance_metric, ' ', num2str(thld),' * STD Surrounding Mean: ', num2str(bad_trials_list(i,1)),'/',num2str(size(EMG_control_all_trials{i,1},1))])
for ii = 1:musc
    % Get muscle trials
    EMG_control_sum = EMG_control_all_trials{i,ii};
    EMG_stim_sum = EMG_stim_all_trials{i,ii};
    % Calc SEM
    SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
    SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));
    % Avg trials
    EMG_control_avg = mean(EMG_control_sum,1);
    EMG_stim_avg = mean(EMG_stim_sum,1);
    % BL subtract
    EMG_control_avg = EMG_control_avg - mean(EMG_control_avg(401:476));
    EMG_stim_avg = EMG_stim_avg - mean(EMG_stim_avg(401:476));

    subplot(1,musc,subplot_pos);
    boundedline(swin-500,EMG_control_avg(swin),SEM_control(swin),'k','alpha')
    hold on
    boundedline(swin-500,EMG_stim_avg(swin),SEM_stim(swin),'c','alpha')
    muscles_names = {'Bicep','Tricep','ECR','PL'};
    animals = {'re7','re8','re9','re10','re12','re14'};
    str = [muscles_names{ii}];%,': Trial Averaged EMG Activty Surrounding Inactivation'];
    title(str)
    ylim([-0.05, 1.6]);
    xlabel('Inactivation at t = 0 ms')
    ylabel('EMG Actvity')
    str1 = ['n = ',int2str(size(EMG_control_sum,1)), ' trials'];
    str2 = ['n = ',int2str(size(EMG_stim_sum,1)), ' trials'];
    %         legend(str1,'control',str2,'stim','Location','southeast')
    hold on

    subplot_pos = subplot_pos + 1;
end

%% Random Sample Algorithm

pseudo_summed_abs_diff_iterations = {};
for iteration = 1:1000
    pseudo_control = {};
    pseudo_stim = {};
    for i = 1:num_animals
        for ii = 1:musc
            % Get random row number
            % Without replacement
            %         pseudo_control_stim_val = randsample(size(EMG_control_all_trials{i,ii},1),size(EMG_control_all_trials{i,ii},1));
            %
            %         pseudo_control_val = pseudo_control_stim_val(1:round(size(pseudo_control_stim_val,1)/2),1);
            %         pseudo_stim_val = pseudo_control_stim_val(round(size(pseudo_control_stim_val,1)/2)+1:size(pseudo_control_stim_val,1),1);

            % With replacement
            pseudo_control_val = randsample(size(EMG_control_all_trials{i,ii},1),round(size(EMG_control_all_trials{i,ii},1)/2),true);
            pseudo_stim_val = randsample(size(EMG_control_all_trials{i,ii},1),round(size(EMG_control_all_trials{i,ii},1)/2),true);

            musc_trials = [];
            for iii = pseudo_control_val
                % Use random row numbers to create separate control and stim trials
                musc_trials = [musc_trials; EMG_control_all_trials{i,ii}(iii,:)];
            end
            pseudo_control{i,ii} = musc_trials;

            musc_trials = [];
            for iiii = pseudo_stim_val
                % Use random row numbers to create separate control and stim trials
                musc_trials = [musc_trials; EMG_control_all_trials{i,ii}(iiii,:)];
            end
            pseudo_stim{i,ii} = musc_trials;
        end
    end

    % Avg trials and change structure of variable
    for i = 1:num_animals
        for ii = 1:musc
            pseudo_control{i,ii} = mean(pseudo_control{i,ii},1);
            pseudo_stim{i,ii} = mean(pseudo_stim{i,ii},1);
        end
    end

    for i = 1:num_animals
        pseudo_control{i,1} = [pseudo_control{i,1};pseudo_control{i,2};pseudo_control{i,3};pseudo_control{i,4}];
        pseudo_stim{i,1} = [pseudo_stim{i,1};pseudo_stim{i,2};pseudo_stim{i,3};pseudo_stim{i,4}];
    end

    for i = 1:3
        pseudo_control(:,2) = [];
        pseudo_stim(:,2) = [];
    end

    % Baseline subtract -20 to 0ms
    pseudo_EMG_control_agg_bl = {};
    pseudo_EMG_stim_agg_bl = {};
    for i = 1:num_animals
        temp_c = [];
        temp_s = [];
        for ii = 1:musc
            temp_c = [temp_c; pseudo_control{i,1}(ii,401:601) - mean(pseudo_control{i,1}(ii,481:501))];
            temp_s = [temp_s; pseudo_stim{i,1}(ii,401:601) - mean(pseudo_stim{i,1}(ii,481:501))];
        end
        pseudo_EMG_control_agg_bl{i,1} = temp_c;
        pseudo_EMG_stim_agg_bl{i,1} = temp_s;
    end

    % Calculate absolute difference per muscle
    pseudo_abs_diff_musc_agg = {};
    for i = 1:num_animals
        pseudo_abs_diff_musc_agg{i, 1} = (abs(pseudo_EMG_control_agg_bl{i, 1} - pseudo_EMG_stim_agg_bl{i, 1}));
    end

        % AVG ACROSS MUSCLES
        % Sum abs diff muscles
        % pseudo_summed_abs_diff = [];
        % for i = 1:num_animals
        %     pseudo_summed_abs_diff = [pseudo_summed_abs_diff; mean(pseudo_abs_diff_musc_agg{i,1})];
        % end

    % AVG ACROSS MICE
    pseudo_summed_abs_diff = {};

    for i = 1:musc
        for j = 1:num_animals
            pseudo_summed_abs_diff{j,i}(1,:) = [pseudo_abs_diff_musc_agg{j,1}(i,:)];
        end
    end

    for i = 1:musc
        for j = 1:num_animals
            pseudo_summed_abs_diff_iterations{j,i}(iteration,:) = pseudo_summed_abs_diff{j,i}(1,:);
        end
    end
end

% Avg iterations per mouse
for i = 1:musc
    for j = 1:num_animals
        pseudo_summed_abs_diff_iterations{j,i} = mean(pseudo_summed_abs_diff_iterations{j,i});
    end
end


disp("Random sample algorithm done!")

%% Average Absolute Difference

% Avg trials and change structure of variable
EMG_control_agg = {};
EMG_stim_agg = {};
for i = 1:num_animals
    for ii = 1:musc
        EMG_control_agg{i,ii} = mean(EMG_control_all_trials{i,ii},1);
        EMG_stim_agg{i,ii} = mean(EMG_stim_all_trials{i,ii},1);
    end
end

for i = 1:num_animals
    EMG_control_agg{i,1} = [EMG_control_agg{i,1};EMG_control_agg{i,2};EMG_control_agg{i,3};EMG_control_agg{i,4}];
    EMG_stim_agg{i,1} = [EMG_stim_agg{i,1};EMG_stim_agg{i,2};EMG_stim_agg{i,3};EMG_stim_agg{i,4}];
end

for i = 1:3
    EMG_control_agg(:,2) = [];
    EMG_stim_agg(:,2) = [];
end

% Baseline subtract -20 to 0ms
EMG_control_agg_bl = {};
EMG_stim_agg_bl = {};
for i = 1:num_animals
    temp_c = [];
    temp_s = [];
    for ii = 1:musc
        temp_c = [temp_c; EMG_control_agg{i,1}(ii,401:601) - mean(EMG_control_agg{i,1}(ii,481:501))];
        temp_s = [temp_s; EMG_stim_agg{i,1}(ii,401:601) - mean(EMG_stim_agg{i,1}(ii,481:501))];
    end
    EMG_control_agg_bl{i,1} = temp_c;
    EMG_stim_agg_bl{i,1} = temp_s;
end

% Calculate absolute difference per muscle
abs_diff_musc_agg = {};
for i = 1:num_animals
    abs_diff_musc_agg{i, 1} = (abs(EMG_control_agg_bl{i, 1} - EMG_stim_agg_bl{i, 1}));
end 

% Sum abs diff muscles
summed_abs_diff = [];
for i = 1:musc
    for j = 1:num_animals
        summed_abs_diff{j,i}(1,:) = [abs_diff_musc_agg{j,1}(i,:)];
    end
end

% Calculate avg absolute difference
avg_diff = [];
for i = 1:num_animals
    % All muscles (averaged then absolute value)
    avg_diff(i,1) = mean(abs_diff_emg_agg(i,111:126)) - mean(abs_diff_emg_agg(i,81:101)); % 25ms window
    avg_diff(i,2) = mean(abs_diff_emg_agg(i,111:151)) - mean(abs_diff_emg_agg(i,81:101)); % 50ms window
    avg_diff(i,3) = mean(abs_diff_emg_agg(i,111:201)) - mean(abs_diff_emg_agg(i,81:101)); % 100ms window

    % All muscles (absolute value then averaged)
    avg_diff(i,4) = mean(summed_abs_diff(i,111:126)) - mean(summed_abs_diff(i,81:101)); % 25ms window
    avg_diff(i,5) = mean(summed_abs_diff(i,111:151)) - mean(summed_abs_diff(i,81:101)); % 50ms window
    avg_diff(i,6) = mean(summed_abs_diff(i,111:201)) - mean(summed_abs_diff(i,81:101)); % 100ms window


    % Individual muscles
    %     avg_diff(i,1) = abs_diff_musc_agg{i,1}(1,111:151) - mean(abs_diff_musc_agg{i,1}(1,81:101)); % Biceps
    %     avg_diff(i,2) = abs_diff_musc_agg{i,1}(2,111:151) - mean(abs_diff_musc_agg{i,1}(2,81:101)); % Triceps
    %     avg_diff(i,3) = abs_diff_musc_agg{i,1}(3,111:151) - mean(abs_diff_musc_agg{i,1}(3,81:101)); % ECR
    %     avg_diff(i,4) = abs_diff_musc_agg{i,1}(4,111:151) - mean(abs_diff_musc_agg{i,1}(4,81:101)); % PL

end

plotmarker = {'o', 's', '^'};
figure; % averaged then absolute value
for i = 1:num_animals/2 % RFA
    scatter([0.75:1:3],avg_diff(i,1:3),75,'r',plotmarker{1,i},'filled');
    hold on
end
for i = (num_animals/2)+1:num_animals % CFA
    scatter([1.25:1:3+0.25],avg_diff(i,1:3),75,'b',plotmarker{1,i-3},'filled');
    hold on
end
for i = 1:3
    plot(i-0.25,mean(avg_diff(1:3,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
    plot(i+0.25,mean(avg_diff(4:6,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
end
x_label = {'10-25ms','10-50ms','10-100ms'};
set(gca,'XTick',1:3,'XTickLabel',x_label, 'Xlimitmethod', 'padded', 'Tickdir', 'out');
ylim([-0.05, 0.3]);
ylabel('Average Difference');
title('Averaged first, then taken abs diff')
% p1 = scatter([0.75:1:3],avg_diff(1,1:3),75,'r', 'filled');
% p2 = scatter([0.75:1:3],avg_diff(2,1:3),75, 'filled');
% p3 = scatter([0.75:1:3],avg_diff(3,1:3),75, 'filled');
% p4 = scatter([1.25:1:3+0.25],avg_diff(1,1:3),75,'b','filled');
% p5 = scatter([1.25:1:3+0.25],avg_diff(2,1:3),75,'filled');
% p6 = scatter([1.25:1:3+0.25],avg_diff(3,1:3),75,'filled');
% legend([p1, p2, p3, p4, p5, p6],animals);
% legend([p1, p4],{'RFA','CFA'},'location','best');
legend(animals,'location','best');

figure; % absolute value then averaged
for i = 1:num_animals/2 % RFA
    scatter([0.75:1:3],avg_diff(i,4:6),75,'r',plotmarker{1,i},'filled');
    hold on
end
for i = (num_animals/2)+1:num_animals % CFA
    scatter([1.25:1:3+0.25],avg_diff(i,4:6),75,'b',plotmarker{1,i-3},'filled');
    hold on
end
for i = 4:6
    plot((i-3)-0.25,mean(avg_diff(1:3,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
    plot((i-3)+0.25,mean(avg_diff(4:6,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
end
x_label = {'10-25ms','10-50ms','10-100ms'};
set(gca,'XTick',1:3,'XTickLabel',x_label, 'Xlimitmethod', 'padded', 'Tickdir', 'out');
ylim([0, 0.3]);
ylabel('Average Difference');
title('Abs diff first, then averaged');
legend(animals,'location','best');

figure;
subplot(1,2,1)
plot([-100:100],summed_abs_diff(1:3,:));
xlabel('Time(ms)');
ylabel('Absolute Difference');
hold on
plot([-100:100],mean(summed_abs_diff(1:3,:),1),'k','LineWidth',2);
xlabel('Time(ms)');
ylabel('Absolute Difference');
ylim([0,0.6]);
hold on
title('RFA Inac');
legend({animals{1:3},'Avg'},'location','best');

subplot(1,2,2)
plot([-100:100],summed_abs_diff(4:6,:));
xlabel('Time(ms)');
ylabel('Absolute Difference');
hold on
plot([-100:100],mean(summed_abs_diff(4:6,:),1),'k','LineWidth',2);
xlabel('Time(ms)');
ylabel('Absolute Difference');
ylim([0,0.6]);
title('CFA Inac');
legend({animals{4:6},'Avg'},'location','best');

% % Organize by muscle
% musc_emg = {};
% for i = 1:musc
%     temp_musc = [];
%     for ii = 1:num_animals
%         temp_musc = [temp_musc; abs_diff_musc_agg{ii,1}(i,:)];
%     end
%     musc_emg{i,1} = temp_musc;
% end

%% Algorithm Part 2

% Subtract algorithm mean with actual data mean
for i = 1:musc
    for j = 1:num_animals
        new_summed_abs_diff{j,i}(1,:) = abs(summed_abs_diff{j,i} - pseudo_summed_abs_diff_iterations{j,i});
        new_summed_abs_diff{j,i}(1,:) = new_summed_abs_diff{j,i}(1,:) - mean(new_summed_abs_diff{j,i}(1,51:101));
    end
end

muscles_names = {'Bicep','Tricep','ECR','PL'};
animals = {'re7','re8','re9','re10','re12','re14'};
figure;
% sgtitle(['Absolute Difference Timeseries-THLD ',num2str(thld),' STD'])
subplot(1,2,1)
plot([-50:100],new_summed_abs_diff(1:3,51:201));
xlabel('Time(ms)');
ylabel('Absolute Difference');
hold on
plot([-50:100],mean(new_summed_abs_diff(1:3,51:201),1),'k','LineWidth',2);
xlabel('Time(ms)');
ylabel('Absolute Difference');
ylim([-0.05, 0.45]);
hold on
title('RFA Inac');
legend({animals{1:3},'Avg'},'location','best');

subplot(1,2,2)
plot([-50:100],new_summed_abs_diff(4:6,51:201));
xlabel('Time(ms)');
ylabel('Absolute Difference');
hold on
plot([-50:100],mean(new_summed_abs_diff(4:6,51:201),1),'k','LineWidth',2);
xlabel('Time(ms)');
ylabel('Absolute Difference');
ylim([-0.05, 0.45]);
title('CFA Inac');
legend({animals{4:6},'Avg'},'location','best');

% Plot avg time series for rfa inac and cfa inac
std_rfa = [];
std_cfa = [];
for i = 1:musc
    std_rfa(i,:) = sqrt(((std([new_summed_abs_diff{1,i}; new_summed_abs_diff{2,i}; new_summed_abs_diff{3,i}],1).^2) .* 3) ./12);
    std_cfa(i,:) = sqrt(((std(new_summed_abs_diff(4:6,:),1).^2) .* 3) ./12);

% std(new_summed_abs_diff(1:3,:),1)
% std(new_summed_abs_diff(4:6,:),1)
SEM_rfa = [];
SEM_cfa = [];
for i = 1:musc
    SEM_rfa(i,:) = (std([new_summed_abs_diff{1,i}; new_summed_abs_diff{2,i}; new_summed_abs_diff{3,i}],1))/sqrt(3);
    SEM_cfa(i,:) = (std([new_summed_abs_diff{4,i}; new_summed_abs_diff{5,i}; new_summed_abs_diff{6,i}],1))/sqrt(3);
end

for i = 1:musc
    figure(1);
    subplot(1,4,i);
    boundedline([-50:100],mean([new_summed_abs_diff{1,i}(51:201); new_summed_abs_diff{2,i}(51:201); new_summed_abs_diff{3,i}(51:201)],1),SEM_rfa(i,51:201),'k','alpha');
    hold on
    boundedline([-50:100],mean([new_summed_abs_diff{4,i}(51:201); new_summed_abs_diff{5,i}(51:201); new_summed_abs_diff{6,i}(51:201)],1),SEM_cfa(i,51:201),'c','alpha');
    % str = animals{1,i};
    % title(['Abssolute Difference Averages-THLD ',num2str(thld),' STD']);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    title([muscles_names{1,i} ' CFA Inactivation']);
    legend('RFA Inac','CFA Inac', 'Location', 'northwest')
    ylim([-0.05,0.45]);
end

for i = 1:4
    figure(1); % line plot
    subplot(1,4,i);
    plot([-50:100],new_summed_abs_diff{i,1}(1,51:201));
    hold on
    plot([-50:100],new_summed_abs_diff{i,2}(1,51:201));
    hold on
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    title([muscles_names{1,i} ' CFA Inactivation']);
    legend('RFA Inac','CFA Inac', 'Location', 'north');
end


% Calculate avg absolute difference
avg_diff = [];
for i = 1:num_animals
   % All muscles (absolute value then averaged)
    avg_diff(i,1) = abs(mean(new_summed_abs_diff(i,111:126)) - mean(new_summed_abs_diff(i,81:101))); % 25ms window
    avg_diff(i,2) = abs(mean(new_summed_abs_diff(i,111:151)) - mean(new_summed_abs_diff(i,81:101))); % 50ms window
    avg_diff(i,3) = abs(mean(new_summed_abs_diff(i,111:201)) - mean(new_summed_abs_diff(i,81:101))); % 100ms window
end

plotmarker = {'o', 'o', 'o'};
figure; % averaged then absolute value
for i = 1:num_animals/2 % RFA
    scatter([0.75:1:3],avg_diff(i,1:3),75,'r',plotmarker{1,i},'filled');
    hold on
end
for i = (num_animals/2)+1:num_animals % CFA
    scatter([1.25:1:3+0.25],avg_diff(i,1:3),75,'b',plotmarker{1,i-3},'filled');
    hold on
end
for i = 1:3
    plot(i-0.25,mean(avg_diff(1:3,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
    plot(i+0.25,mean(avg_diff(4:6,i)),'k_','MarkerSize',12,'LineWidth',2);
    hold on
end
x_label = {'10-25ms','10-50ms','10-100ms'};
set(gca,'XTick',1:3,'XTickLabel',x_label, 'Xlimitmethod', 'padded', 'Tickdir', 'out');
ylim([0, 0.25]);
ylabel('Average Difference');
% title(['Average Absolute Difference-THLD ',num2str(thld),' STD'])
% legend(animals,'location','best');

%% Significance Calculations
rfa_abs_diff_mean = mean(new_summed_abs_diff(1:3,51:201),1);
cfa_abs_diff_mean = mean(new_summed_abs_diff(4:6,51:201),1);

% mean abs diff from -51-0ms
rfa_abs_diff_zero = mean(rfa_abs_diff_mean(1,1:51));
cfa_abs_diff_zero = mean(cfa_abs_diff_mean(1,1:51));

% abs diff data set from 0-25/50/100ms
rfa_abs_diff_25 = rfa_abs_diff_mean(1,51:76);
rfa_abs_diff_50 = rfa_abs_diff_mean(1,51:101);
rfa_abs_diff_100 = rfa_abs_diff_mean(1,51:151);

cfa_abs_diff_25 = cfa_abs_diff_mean(1,51:76);
cfa_abs_diff_50 = cfa_abs_diff_mean(1,51:101);
cfa_abs_diff_100 = cfa_abs_diff_mean(1,51:151);

% paired t-tests 
[rfa_sig_25, rfa_sig_25_p] = ttest(rfa_abs_diff_25,rfa_abs_diff_zero, 'Tail', 'right')
[rfa_sig_50, rfa_sig_50_p] = ttest(rfa_abs_diff_50,rfa_abs_diff_zero, 'Tail', 'right')
[rfa_sig_100, rfa_sig_100_p] = ttest(rfa_abs_diff_100,rfa_abs_diff_zero, 'Tail', 'right')

[cfa_sig_25, cfa_sig_25_p] = ttest(cfa_abs_diff_25,cfa_abs_diff_zero, 'Tail', 'right')
[cfa_sig_50, cfa_sig_50_p] = ttest(cfa_abs_diff_50,cfa_abs_diff_zero, 'Tail', 'right')
[cfa_sig_100, cfa_sig_100_p] = ttest(cfa_abs_diff_100,cfa_abs_diff_zero, 'Tail', 'right')

sigvalues = [rfa_sig_25_p, cfa_sig_25_p; rfa_sig_50_p, cfa_sig_50_p; rfa_sig_100_p, cfa_sig_100_p];


%% Muscle Effects Per Mouse

figure; % Trial Averaged EMG Activty Surrounding Laser Inactivation
sgtitle('Biceps EMG Activity'); 
for i = 1:num_animals
    subplot(2,num_animals,i)
    boundedline([-100:1:100],EMG_control_agg{i,1}(1,:),SEM_control_agg{i,1}(1,:),'k','alpha');
    hold on
    boundedline([-100:1:100],EMG_stim_agg{i,1}(1,:),SEM_stim_agg{i,1}(1,:),'c','alpha');
    str = animals{1,i};
    title(str);
%     xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([-0.5,3.5]);
    hold on
end

for i = 1+num_animals:num_animals*2
    subplot(2,num_animals,i)
    plot([-100:1:100],abs_diff_musc_agg{i-num_animals,1}(1,:),'Color','k');
    hold on
%     str = animals{1,i};
%     title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    ylim([-0.1,0.6]);
    ii = ii + 1;
    hold on
end

figure; % Trial Averaged EMG Activty Surrounding Laser Inactivation
sgtitle('Triceps EMG Activity'); 
for i = 1:num_animals
    subplot(2,num_animals,i)
    boundedline([-100:1:100],EMG_control_agg{i,1}(2,:),SEM_control_agg{i,1}(2,:),'k','alpha');
    hold on
    boundedline([-100:1:100],EMG_stim_agg{i,1}(2,:),SEM_stim_agg{i,1}(2,:),'c','alpha');
    str = animals{1,i};
    title(str);
%     xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([-0.5,3.5]);
    hold on
end

for i = 1+num_animals:num_animals*2
    subplot(2,num_animals,i)
    plot([-100:1:100],abs_diff_musc_agg{i-num_animals,1}(2,:),'Color','k');
    hold on
%     str = animals{1,i};
%     title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    ylim([-0.1,0.6]);
    ii = ii + 1;
    hold on
end

figure; % Trial Averaged EMG Activty Surrounding Laser Inactivation
sgtitle('ECR EMG Activity'); 
for i = 1:num_animals
    subplot(2,num_animals,i)
    boundedline([-100:1:100],EMG_control_agg{i,1}(3,:),SEM_control_agg{i,1}(3,:),'k','alpha');
    hold on
    boundedline([-100:1:100],EMG_stim_agg{i,1}(3,:),SEM_stim_agg{i,1}(3,:),'c','alpha');
    str = animals{1,i};
    title(str);
%     xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([-0.5,3.5]);
    hold on
end

for i = 1+num_animals:num_animals*2
    subplot(2,num_animals,i)
    plot([-100:1:100],abs_diff_musc_agg{i-num_animals,1}(3,:),'Color','k');
    hold on
%     str = animals{1,i};
%     title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    ylim([-0.1,0.6]);
    ii = ii + 1;
    hold on
end

figure; % Trial Averaged EMG Activty Surrounding Laser Inactivation
sgtitle('PL EMG Activity'); 
for i = 1:num_animals
    subplot(2,num_animals,i)
    boundedline([-100:1:100],EMG_control_agg{i,1}(4,:),SEM_control_agg{i,1}(4,:),'k','alpha');
    hold on
    boundedline([-100:1:100],EMG_stim_agg{i,1}(4,:),SEM_stim_agg{i,1}(4,:),'c','alpha');
    str = animals{1,i};
    title(str);
%     xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([-0.5,3.5]);
    hold on
end

for i = 1+num_animals:num_animals*2
    subplot(2,num_animals,i)
    plot([-100:1:100],abs_diff_musc_agg{i-num_animals,1}(4,:),'Color','k');
    hold on
%     str = animals{1,i};
%     title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    ylim([-0.1,0.6]);
    ii = ii + 1;
    hold on
end

figure; % Trial Averaged EMG Activty Surrounding Laser Inactivation
sgtitle('Summed Muscle EMG Activity'); 
for i = 1:num_animals
    subplot(2,num_animals,i)
    boundedline([-100:1:100],EMG_control_agg_musc_avg{i,1}(1,:),SEM_control_agg_musc_avg{i,1}(1,:),'k','alpha');
    hold on
    boundedline([-100:1:100],EMG_stim_agg_musc_avg{i,1}(1,:),SEM_stim_agg_musc_avg{i,1}(1,:),'c','alpha');
    str = animals{1,i};
    title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([-0.5,5.5]);
    hold on
end

ii = 1;
for i = 1+num_animals:num_animals*2
    subplot(2,num_animals,i)
    plot([-100:1:100],abs_diff_emg_agg(ii,:),'Color','k');
    hold on
%     str = animals{1,i};
%     title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('Absolute Difference');
    ylim([-0.1,0.6]);
    ii = ii + 1;
    hold on
end

%%

RFA_control_avg = mean(cell2mat(EMG_control_agg_musc_avg(1:3,1)),1);
CFA_control_avg = mean(cell2mat(EMG_control_agg_musc_avg(3:6,1)),1);
RFA_stim_avg = mean(cell2mat(EMG_stim_agg_musc_avg(1:3,1)),1);
CFA_stim_avg = mean(cell2mat(EMG_stim_agg_musc_avg(3:6,1)),1);
RFA_SEM_control = mean(cell2mat(SEM_control_agg_musc_avg(1:3,1)),1);
CFA_SEM_control = mean(cell2mat(SEM_control_agg_musc_avg(3:6,1)),1);
RFA_SEM_stim = mean(cell2mat(SEM_stim_agg_musc_avg(1:3,1)),1);
CFA_SEM_stim = mean(cell2mat(SEM_stim_agg_musc_avg(3:6,1)),1);

figure();
subplot(2,1,2) % RFA
boundedline([-100:1:100],RFA_control_avg,RFA_SEM_control,'k','alpha');
hold on
boundedline([-100:1:100],RFA_stim_avg,RFA_SEM_stim,'c','alpha');
% str = animals{1,i};
title('RFA Inactivation');
xlabel('Inactivation at t = 0 ms');
ylabel('EMG Actvity');
ylim([-1,3]);
hold on

subplot(2,1,1) % CFA
boundedline([-100:1:100],CFA_control_avg,CFA_SEM_control,'k','alpha');
hold on
boundedline([-100:1:100],CFA_stim_avg,CFA_SEM_stim,'c','alpha');
% str = animals{1,i};
title('CFA Inactivation');
xlabel('Inactivation at t = 0 ms');
ylabel('EMG Actvity');
ylim([-1,3]);
hold on

%%

RFA_abs_diff = (abs(RFA_control_avg - RFA_stim_avg));
CFA_abs_diff = (abs(CFA_control_avg - CFA_stim_avg));
RFA_SEM_abs_diff = (abs(RFA_SEM_control - RFA_SEM_stim));
CFA_SEM_abs_diff = (abs(CFA_SEM_control - CFA_SEM_stim));


figure();
boundedline([-100:1:100],RFA_abs_diff, RFA_SEM_abs_diff,'k','alpha');
% plot([-100:1:100], RFA_abs_diff);
hold on
boundedline([-100:1:100],CFA_abs_diff, CFA_SEM_abs_diff,'c','alpha');
% plot([-100:1:100], CFA_abs_diff);
p1 = boundedline([-100:1:100],RFA_abs_diff, RFA_SEM_abs_diff,'k','alpha');
p2 = boundedline([-100:1:100],CFA_abs_diff, CFA_SEM_abs_diff,'c','alpha');
legend([p1, p2],{'RFA Inac','CFA Inac'}, 'location','best');

%% plot bins
% 
% for ii = 1:num_neurons
%     figure;
%     bin_stim_trials = squeeze(bin_stim(:,400:600,ii));
%     subplot(2,1,1);
%     imagesc(bin_stim_trials)
%     title('stim')
%     bin_control_trials = squeeze(bin_control(:,400:600,ii));
%     subplot(2,1,2);
%     imagesc(bin_control_trials)
%     title('control')
%     str = num2str(depths(ii));
%     sgtitle(str)
% end

%%

figure(1);
for i = 1:size(success_failed_reaches,2)
    subplot(4,6,i)
    bar([1:4],[success_failed_reaches{3,i};success_failed_reaches{4,i}])
    title([success_failed_reaches{1,i}, '-', success_failed_reaches{2,i}])
    xlabel('Spout')
    ylim([0,100])
end
legend('success','failed')
sgtitle('Success/Failure Reaches')

figure(2);
for i = 1:size(success_failed_reaches,2)
    subplot(4,6,i)
    bar((success_failed_reaches{3,i} ./ (success_failed_reaches{3,i} + success_failed_reaches{4,i})))
    hold on
    title([success_failed_reaches{1,i}, '-', success_failed_reaches{2,i}])
    xlabel('Spout')
    ylim([0,1])
end
sgtitle('Success Reaches Percentage')
    
