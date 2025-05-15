%% INPUT 'RFA', 'CFA', or 'ALL'
animal_subset = 'NEW'

%%
addpath(genpath('Z:\Scripts'));
if strcmp(animal_subset,'RFA')
    animals = {'re7','re8','re9'};
elseif strcmp(animal_subset,'CFA')
    animals = {'re10','re12'};
elseif strcmp(animal_subset,'ALL')
    animals = {'re7','re8','re9','re10','re12'};
elseif strcmp(animal_subset,'NEW')
    animals = {'re7','re9'};
else
    disp('Error! Wrong aniamls!')
end
num_animals = length(animals);
abs_diff_emg_agg = [];

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
        str = ['Z:\reaching\',animal,'_',date,'_g0\preprocess\EMGoneSided.mat'];
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
%         figure;
%         boundedline(time(swin),EMG_control_sub(swin),SEM_control(swin),'k','alpha')
%         hold on
%         boundedline(time(swin),EMG_stim_sub(swin),SEM_stim(swin),'c','alpha')
%         str = [animal,' - ',date,': Trial Averaged EMG Activty Surrounding Laser Inactivation'];
%         title(str)
%         xlabel('Inactivation at t = 0 ms')
%         ylabel('EMG Actvity')
%         str1 = ['n = ',int2str(length(control_times)), ' trials'];
%         str2 = ['n = ',int2str(length(stim_times)), ' trials'];
%         legend(str1,'control',str2,'stim','Location','southeast')
% 
%         figure;
%         boundedline(time(swin),fr_control_sub(swin),fr_SEM_control(swin),'k','alpha')
%         hold on
%         boundedline(time(swin),fr_stim_sub(swin),fr_SEM_stim(swin),'c','alpha')
%         str = [animal,' - ',date,': Trial Averaged CFA Activty Surrounding Laser Inactivation'];
%         title(str)
%         xlabel('Inactivation at t = 0 ms')
%         ylabel('Neural Actvity')
%         str1 = ['n = ',int2str(length(control_times)), ' trials'];
%         str2 = ['n = ',int2str(length(stim_times)), ' trials'];
%         legend(str1,'control',str2,'stim','Location','southeast')

   
    end


%% Individual and Aggregate Muscle Activity
EMG_control_agg_mouse = [];
EMG_stim_agg_mouse = [];
for musc = 1:num_muscles
    EMG_control_sum = EMG_all_control(:,:,musc);
    EMG_stim_sum = EMG_all_stim(:,:,musc);
    
    SEM_control = std(EMG_control_sum,1)/sqrt(size(EMG_control_sum,1));
    SEM_stim = std(EMG_stim_sum,1)/sqrt(size(EMG_stim_sum,1));
    
    EMG_control_trial_avg = mean(EMG_control_sum);
    EMG_stim_trial_avg = mean(EMG_stim_sum);
    EMG_control_bl_mean = mean(EMG_control_trial_avg(1:300));
    EMG_stim_bl_mean = mean(EMG_stim_trial_avg(1:300));
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

%     %figure;
%     subplot(1,4,musc);
%     boundedline(time(swin),EMG_control_sub(swin),SEM_control(swin),'k','alpha')
%     hold on
%     boundedline(time(swin),EMG_stim_sub(swin),SEM_stim(swin),'c','alpha')
%     muscles_names = {'Bicep','Tricep','ECR','PL'};
%     str = [animal,' - ',muscles_names{musc}];%,': Trial Averaged EMG Activty Surrounding Inactivation'];
%     %title(str)
%     ylim([-0.1 2]);
%     %xlabel('Inactivation at t = 0 ms')
%     %ylabel('EMG Actvity')
%     str1 = ['n = ',int2str(size(EMG_control_sum,1)), ' trials'];
%     str2 = ['n = ',int2str(size(EMG_stim_sum,1)), ' trials'];
%     %legend(str1,'control',str2,'stim','Location','southeast')
end

% Aggregate muscle activity from all animals and muscles
EMG_control_agg{animal_index,:} = EMG_control_agg_mouse;
EMG_stim_agg{animal_index,:} = EMG_stim_agg_mouse;
SEM_control_agg{animal_index,:} = SEM_control_agg_mouse;
SEM_stim_agg{animal_index,:} = SEM_stim_agg_mouse;
end

%% Avg EMG across mice broke out by muscle

% Avg control and stim
musc_traces_control = {};
musc_traces_stim = {};
for i = 1:musc
    temp_c = [];
    temp_s = [];
    for ii = 1:num_animals
        temp_c = [temp_c; EMG_control_agg{ii,1}(i,:)];
        temp_s = [temp_s; EMG_stim_agg{ii,1}(i,:)];
    end
    musc_traces_control{i,1} = temp_c;
    musc_traces_stim{i,1} = temp_s;
end
musc_traces_avg_control = [];
musc_traces_avg_stim = [];
for i = 1:musc
    musc_traces_avg_control(i,:) = mean(musc_traces_control{i,1},1);
    musc_traces_avg_stim(i,:) = mean(musc_traces_stim{i,1},1);
end

% SEM calc
SEM_musc_traces_control = [];
SEM_musc_traces_stim = [];
for i = 1:musc
    SEM_musc_traces_control(i,:) = std(musc_traces_control{i,1},1)/sqrt(size(musc_traces_control{i,1},1));
    SEM_musc_traces_stim(i,:) = std(musc_traces_stim{i,1},1)/sqrt(size(musc_traces_stim{i,1},1));

end

% Plot control
musc_label = {'Biceps','Triceps','ECR','PL'};
figure();
sgtitle('EMG Control Trials Averages'); 
for i = 1:musc
    subplot(musc,1,i)
    boundedline([-100:1:100],musc_traces_avg_control(i,:),SEM_musc_traces_control(i,:),'k','alpha');
    str = musc_label{1,i};
    title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([0,2]);
    hold on
end

% Plot stim
figure();
sgtitle('EMG Stim Trials Averages'); 
for i = 1:musc
    subplot(musc,1,i)
    boundedline([-100:1:100],musc_traces_avg_stim(i,:),SEM_musc_traces_stim(i,:),'k','alpha');
    str = musc_label{1,i};
    title(str);
    xlabel('Inactivation at t = 0 ms');
    ylabel('EMG Actvity');
    ylim([0,2]);
    hold on
end




























