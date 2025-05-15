 %% Trial Muscle Averages
 
 %%%%%% DESCRIPTION %%%%%%
 % This script takes the output from proprocessing_with_acg.m and creates
 % trial bounds and pyramidal neuron data structures.
 %%%%%%%%%%%%%%%%%%%%%%%%%
 
%% Add Paths
% Add paths to needed functions in this script
% (These paths need to be edited to the proper location of these functions)
opengl software
home_path = 'Z:\Scripts';
addpath(home_path);
addpath('Z:\neuropixel');
neur_paths = genpath('Z:\Scripts');
addpath(neur_paths);

%% Global Variables

% Recording Identifiers
dates = ['03272020';'03302020';'03312020';'04012020';'04022020';...
    '04032020';'09222020';'09232020';'09302020';'10012020';'10022020';...
    '10052020';'10062020';'10242020';'10262020';'10272020';'10282020';...
    '11102020';'11122020';'11152020';'11172020'];
dates = ['03302020'];

animals = {'a048';'a048';'a048';'a048';'a048';'a048';'a050';'a050';...
    'ss2';'ss2';'ss2';'ss2';'ss2';'MA1';'MA1';'MA1';'MA1';...
    'a051';'a051';'MA2';'MA2'};
animals = {'a048'};

all_RFA_cutoffs = {'0327',4000;'0330',4000;'0331',4000;'0401',4000;'0402',4000;...
                   '0403',4000;'0922',2100;'0923',2100;'0930',2100;'1001',2100;...
                   '1002',3000;'1005',2200;'1006',2200;'1024',2000;'1026',2500;...
                   '1027',2100;'1028',2100;'1110',2500;'1112',2500;'1115',2100;...
                   '1117',2100;};

% Function Variables
eliminate_bad = 1; %0 - do not eliminate pathological neurons in RFA
                   %1 - eliminate pathological neurons in RFA
should_save = 0; %0 - do not save ...
                 %1 - save ...
plot_all_EMG = 0;
plot_FR_trialavg = 0;
                  
% Digital channel numbers
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

% Addtional Variables
num_spouts = 4;
lift_window = [-200 800]; 
kilosort_dir = 'kilosort';
preprocess_dir = 'preprocess_with_acg';



%%

for date_number = 1:size(dates,1)  
% Set max index for RFA neuron

    date = dates(date_number,:);
    animal = animals{date_number};
    fprintf('recording %s, animal %s\n',date,animal);

    % Set maximum index for removal of pathological RFA neurons assigned to
    % depths outside of tissue.
    total_num_sessions = size(all_RFA_cutoffs,1);
    for ii = 1:total_num_sessions
        if all_RFA_cutoffs{ii,1}==date(1:4)
            session_index = ii;
            break
        end
        if ii == total_num_sessions
            disp('Error! Date provided is incorrect')
        end
    end
    max_good_index = all_RFA_cutoffs{session_index,2};

%% Load Data

    str = sprintf('Z:\\akiko\\%s_%s_g0\\%s\\',date,animal,preprocess_dir);
    cd(str);
    files = dir([str '*.mat']);

    % Loads .mat files from selected directory
    for i=1:length(files)
        load(files(i).name,'-mat');
    end
    cd(home_path);

    % Uses EMG file to determine number of muscles recorded from
    num_muscles = size(EMG,1);
    if (num_muscles == 4)
        muscles = ["Bicep","Tricep","ECR","PL"];
    elseif (num_muscles == 6)
         muscles = ["Deltoid","Pectoralis","Bicep","Tricep","ECR","PL"]; 
    end

    disp('done loading data!');

%% Remove ISI violations

    % Important variables: duration of the time series (in ms) and the
    % num_probes used in this data set.
    duration = size(EMG,2);
    num_probes = length(events);
    num_all_good(1) = length(unique(events{1}(:,1))); 
    num_all_good(2) = length(unique(events{2}(:,1)));
    events = removeFlaggedUnits(events,unit_flags);
    for ii = 1:num_probes
        fprintf('Probe %d: %d ISI violations out of %d total good units\n',...
            ii-1,length(unit_flags{ii}),num_all_good(ii));
    end
    
%% Generate reward and units cells

    % Create reward cell containing the times when the reward sound is
    % detected (this is for both successful and unsuccesful trials).
    rw = get_reward_times(digital, [reward1 reward2 reward3 reward4]);

    % correct for the reward sound by subtracting 100
    rw1 = rw{reward1}-100;
    rw2 = rw{reward2}-100;
    rw3 = rw{reward3}-100;
    rw4 = rw{reward4}-100;
    rw_corr = {rw1,rw2,rw3,rw4};

    % Make empty cells for the train, FR, and units.
    trains = cell(num_probes,1);
    units = cell(num_probes,1);
    FR = cell(num_probes,1);
    num_units = zeros(num_probes,1);

    % Iterate across the number of probes to generate the trains and units
    for i=1:num_probes

        % kilosort starts indexing at 0 but lets just add 1 to the index to
        % make things easier
        units{i} = unique(events{i}(:,1));
        if (min(units{i})==0)
            disp('indexing starts at 0, fixing that now...');
            events{i}(:,1) = events{i}(:,1)+1;
            units{i} = units{i}+1;
        end

        % events_to_train() takes the events matrix and outputs a 
        % (num_probe x 1) cell in which each cell contains a 
        % (length(units) x 1) matrix of all the spike times for a given unit
        trains{i} = events_to_train(events{i});

        % train_to_FR() takes the trains cell and outputs a
        % (num_probe x 1) cell in which each cell contains a 
        % (length(units) x duration) matrix where each row is a smoothed 
        % time series of the firing rate of that unit.
        FR{i} = train_to_FR(trains{i},units{i});

        % Sometimes FR has an extra time point at the end. We remove it so that
        % both EMG and FR have a length of duration.
        d = size(FR{i},2) - size(EMG,2);   
        FR{i} = FR{i}(:,1:end-d);

        num_units(i) = length(units{i});  
    end

%% Get CFA/RFA indices

    % The depth of 1500 was determined from previous experiments. Events that
    % occur at depths superior to the threshold are assumed to represent cortex
    % units. events_to_cor isolates the units that are CFA (probe 1) and RFA
    % (probe 2). Note that the contents of CFA_RFA_units are the unit names
    % themselves (i.e. NOT the index of the unit name).

    if eliminate_bad == 1
        bad_events_inds = find(events{2}(:,3)>max_good_index);
        bad_unit_inds = unique(events{2}(bad_events_inds,1));
        events{2}(bad_events_inds,:) = [];
        for ii = 1:length(bad_unit_inds)
            units{2}(find(units{2}==bad_unit_inds(ii))) = [];
        end
        str = sprintf('CFA: 0 units Out of Tissue');
        disp(str)
        str = sprintf('RFA: %d units Out of Tissue',length(bad_unit_inds));
        disp(str)
    else
        str = sprintf('CFA: 0 units Out of Tissue');
        disp(str)
        str = sprintf('RFA: 0 units Out of Tissue');
        disp(str)
    end
    nogo = 0;
    if length(animal) == 4
        if animal(1) == 'a' && animal(4) == '8'
            limit_low = 1500;
            limit_high = 0;
            [~,CFA_RFA_units] = ...
                events_to_cor_old(num_probes,events,limit_low,limit_high);
            CFA_lower_depth = max(events{1}(:,3))-limit_low;
            RFA_lower_depth = max(events{2}(:,3))-limit_low;
            nogo = 1;
        end
    end
    if nogo == 0
        limit_low = 1500;
        [CFA_train,RFA_train,CFA_lower_depth,RFA_lower_depth] = ...
            events_to_cor(events,limit_low,max_good_index);
        CFA_RFA_trains = {CFA_train,RFA_train};
        CFA_RFA_units = cell(num_probes,1);
        for ii = 1:num_probes
            for jj = 1:length(CFA_RFA_trains{ii})
                if ~isempty(CFA_RFA_trains{ii}{jj})
                    CFA_RFA_units{ii} = [CFA_RFA_units{ii};jj];
                end
            end
        end
    end
    CFA_units = CFA_RFA_units{1};
    RFA_units = CFA_RFA_units{2};

    all_depths = cell(num_probes,1);
    for ii = 1:num_probes
        all_depths{ii} = zeros(length(units{ii}),1);
        for jj = 1:length(units{ii})
            unit_ind = find(events{ii}(:,1)==units{ii}(jj),1);
            all_depths{ii}(jj) = events{ii}(unit_ind,3);
        end
        figure;
        histogram(all_depths{ii},40);
        str = sprintf('Probe %d: Neuron Depths (all)',ii-1);
        title(str);
        xlabel('Depth (from tip of probe)')
        if ii == 1
            xline(CFA_lower_depth,'r','CFA cutoff','LineWidth',2.0)
        elseif ii == 2
            xline(RFA_lower_depth,'r','RFA cutoff','LineWidth',2.0)
        end
    end

    % Display how many units on each probe are being analyzed
    for ii=1:num_probes
        str = sprintf('We have %d good units on probe %d',length(units{ii}),ii-1);
        disp(str)
    end

    %% Getting reach data
    % Creates cell touch_channels that contains digital data for all of the
    % beam breaks at the spout.
    touch_channel1 = digital(touch_sensor1,:);
    touch_channel2 = digital(touch_sensor2,:);
    touch_channel3 = digital(touch_sensor3,:);
    touch_channel4 = digital(touch_sensor4,:);
    touch_channels = [touch_channel1; touch_channel2; ...
        touch_channel3; touch_channel4];

    % find_success() creates a series of time bounds for every successful
    % trial. A succesful trial is defined as a touch (beam break) to the
    % correct spout before suction, providing that the correct spout is touched
    % within the allowable_error after an incorrect spout. The time bounds are
    % given by [reward time, touch time].
    allowable_error = 50;
    success_bounds = find_success(num_spouts,rw_corr,touch_channels,...
        digital(suction,:),allowable_error);

    all_success_nospout = [];
    for ii = 1:num_spouts
        all_success_nospout = [all_success_nospout; success_bounds{ii}];
    end

    % inact_muscles are the muscles that allow us to observe reasonably flat
    % baselines of EMG activty (the shoulder muscles are often active in the
    % quiescence period
    inact_muscles = muscles~='Deltoid' & muscles~='Pectoralis';

    % reach_bounds() creates a series of time bounds [start of reach, touch]
    % for each successful trial.
    [reach_bounds, EMG_reach, rel_reach_times, up_thresh_times, peak_mag,...
        peak_times, ramp, trial_std] = get_reach_bounds(date,num_spouts,...
        EMG(inact_muscles,:),success_bounds);

    % create max indices
    spout_max_index = cell(1,num_spouts);
    running = 0;
    for ii = 1:num_spouts
        spout_max_index{ii} = running + length(reach_bounds{ii}(:,1));
        running = spout_max_index{ii};
    end

    % Flatten reaches to remove spout dependence.
    all_reach_nospout = [];
    for ii = 1:num_spouts
        all_reach_nospout = [all_reach_nospout; reach_bounds{ii}];
    end

    % Flatten ramp times to remove spout dependence.
    all_ramp_nospout = [];
    for ii = 1:num_spouts
        all_ramp_nospout = [all_ramp_nospout, ramp{ii}];
    end


%% Make Histogram of Touch Latency (touch - reach)

    % Determine the latency of the touch compared to the reach for each
    % successful trial. The latencies are separated in terms of spout in order
    % to determine if there is spout dependence on the latency.
    latency = cell(num_spouts,1);
    for ii = 1:num_spouts
        spout_trials = length(reach_bounds{ii}(:,1));
        latency{ii} = zeros(spout_trials,1);
        for jj = 1:spout_trials
            latency{ii}(jj) = reach_bounds{ii}(jj,2)-reach_bounds{ii}(jj,1);         
        end
    end

    % Generates and displays the latency histogram for each spout and then
    % combines all latencies into a single matrix.
    all_latency = [];
    for ii = 1:num_spouts
        all_latency = [all_latency,latency{ii}'];
    end

    % Calculate and display the mean and std for all latencies. 
    latency_mean = mean(all_latency);
    latency_std = std(all_latency);
% 
%     figure;
%     histogram(all_latency,30);
%     xline(latency_mean+(3*latency_std));
%     title('Touch Latency (touch - reach)');

%% Make Histogram of Reach Latency (reach - reward)

    % Determine the latency of the reach compared to the reward for each
    % successful trial. The latencies are separated in terms of spout in order
    % to determine if there is spout dependence on the latency.
    reward_latency = cell(num_spouts,1);
    for ii = 1:num_spouts
        spout_trials = length(success_bounds{ii}(:,1));
        reward_latency{ii} = zeros(spout_trials,1);
        for jj = 1:spout_trials
            reward_latency{ii}(jj) = reach_bounds{ii}(jj,1)-...
                success_bounds{ii}(jj,1);
        end
    end

    % Generates and displays the latency histogram for each spout and then
    % combines all latencies into a single matrix.
    all_reward_latency = [];
    for ii = 1:num_spouts
         all_reward_latency = [all_reward_latency,reward_latency{ii}'];
    end

    % Calculate and display the mean and std for all latencies. 
    reward_latency_mean = mean(all_reward_latency);
    reward_latency_std = std(all_reward_latency);

%     figure;
%     histogram(all_reward_latency,30);
%     title('Reach Latency (reach - reward)');


%% Flag Indices

    %Display Total Number of trials
    str = sprintf('Total Number of Successful Reaches: %d',length(all_latency));
    disp(str)

    %Flag those with >mean+3*std ramp times.
    flag_ind_ramp = find(all_ramp_nospout>...
        mean(all_ramp_nospout+3*std(all_ramp_nospout)));
    str = sprintf('%d trials flagged for ramp time > mean+3*std.'...
        ,length(flag_ind_ramp));
    disp(str)

    figure;
    histogram(all_ramp_nospout,30);
    xline(mean(all_ramp_nospout+3*std(all_ramp_nospout)));
    title('Ramp Times (reach to touch)');

    %Flag those with reach time <450ms
    all_rel_reach_times = (all_reach_nospout(:,1) - all_success_nospout(:,1))';
    flag_ind_early = find(all_rel_reach_times<-50);
    str = sprintf('%d trials flagged for reach time prior to -50ms.'...
        ,length(flag_ind_early));
    disp(str)

    % figure;
    % histogram(all_rel_reach_times,30);
    % xline(-50);
    % title('Reach Latency (reach - reward)');

    % Flag indices corresponding to 'zero' latencies
    flag_ind_zero = find(all_latency==0);
    str = sprintf('%d trials flagged for lack of threshold crossing.'...
        ,length(flag_ind_zero));
    disp(str)

    % Flag indices with > mean+3*std reach/touch latencies
    flag_ind_long = find(all_latency>latency_mean+(3*latency_std));
    str = sprintf('%d trials flagged for long reach/touch latency.'...
        ,length(flag_ind_long));
    disp(str)

    % Flag indicies with > mean+3*std of standard deviations of intial 
    % baselines
    flag_ind_noisy = find(trial_std>mean(trial_std)+3*std(trial_std));
    str = sprintf('%d trials flagged for noisy baselines.'...
        ,length(flag_ind_noisy));
    disp(str)

    figure;
    histogram(trial_std,30);
    xline(mean(trial_std)+3*std(trial_std));
    title('Individual Trial Baseline St. Dev.');

    % Concatenate all flags into one matrix
    all_flags = [flag_ind_early,flag_ind_zero,flag_ind_ramp,...
        flag_ind_long,flag_ind_noisy];
    all_flags = unique(all_flags);

    % Display the number of unique trials flagged to removed
    str = sprintf('There were %d unique flags',length(all_flags));
    disp(str)

    % find number of flags for each spout
    flags_byspout = cell(1,num_spouts);
    for ii = 1:num_spouts
        if ii == 1
            flags_byspout{ii} = length(find(all_flags<=spout_max_index{ii}));
        else
            flags_byspout{ii} = length(find(all_flags>spout_max_index{ii-1} &...
                all_flags<=spout_max_index{ii}));
        end
    end

    % recalculate the max index
    running = 0;
    spout_max_index_edit = cell(1,num_spouts);
    for ii = 1:num_spouts
        spout_max_index_edit{ii} = spout_max_index{ii}-flags_byspout{ii}-running;
        running = running + flags_byspout{ii};
    end

%% Make Edited Reach Bounds and Latency Arrays and Histograms

    % Removed flagged trials from all relevent matrices.
    EMG_reach_edit = EMG_reach;
    EMG_reach_edit(all_flags) = [];

    reach_bounds_edit = all_reach_nospout;
    reach_bounds_edit(all_flags,:) = [];

    up_thresh_times_edit = up_thresh_times;
    up_thresh_times_edit(all_flags) = [];

    rel_reach_times_edit = rel_reach_times;
    rel_reach_times_edit(all_flags) = [];

    peak_times_edit = peak_times;
    peak_times_edit(all_flags) = [];

    peak_mag_edit = peak_mag;
    peak_mag_edit(all_flags) = [];

    all_latency_edit = all_latency;
    all_latency_edit(all_flags) = [];

    all_reward_latency_edit = all_reward_latency;
    all_reward_latency_edit(all_flags) = [];

    all_ramp_edit = all_ramp_nospout;
    all_ramp_edit(all_flags) = [];

    trial_std_edit = trial_std;
    trial_std_edit(all_flags) = [];
    
%% Plot EMG traces (red if flagged, blue if otherwise)
    if plot_all_EMG == 1
        for ii = 1:length(rel_reach_times_edit)
            figure;     
            plot(EMG_reach_edit{ii}); 
            hold on
            plot(rel_reach_times_edit(ii),...
                EMG_reach_edit{ii}(rel_reach_times_edit(ii)),'d');
        end

        for ii = 1:length(rel_reach_times)
            if ismember(ii,all_flags)
                figure;     
                plot(EMG_reach{ii},'r'); 
                hold on
                plot(rel_reach_times(ii),EMG_reach{ii}(rel_reach_times(ii)),'d');
            end
        end
    end
%%
    % Display the number of reach trials used for analysis
    num_good_trials = length(rel_reach_times_edit);
    str = sprintf('Number of Used Trials = %d',num_good_trials);
    disp(str)

    % Display histogram of relative reach times (should range between -50ms to
    % ~600ms).
    figure;
    histogram(rel_reach_times_edit-500,30);
    react_mean = mean(rel_reach_times_edit-500);
    str = sprintf('Reach Latency - Edited (reach - reward): Mean = %5.1f',react_mean);
    title(str);
    xlabel('Miliseconds After Reward');
    ylabel('Number of Trials');

%% Make Lift Bounds in order to Trim FR matrix

    % We create windows -150ms to +50ms surrounding the reach startclear
     
    lift_duration = lift_window(2)-lift_window(1)+1;
    lift_bounds = zeros(length(reach_bounds_edit(:,1)),2);

    % lift_bounds is a (number of successful trials for all spouts x 2) matrix
    % that contains all of the lift windows for all accepted trials.
    for ii = 1:length(reach_bounds_edit(:,1))
        lift_bounds(ii,:) = [reach_bounds_edit(ii,1)+lift_window(1),...
           reach_bounds_edit(ii,1)+lift_window(2)];
    end

%% Take Trial Averages of Firing Rates During Liftbounds

    % Trim the FR matrix to only the windows around the reach initiation
    % Arranged in a (num_units x lift_duration x total_num_trials) matrix 
    % containing the FR time series for each unit for every trial.
    all_trials = cell(num_probes,1);
    for ii = 1:num_probes
        for jj = 1:length(units{ii}(:,1))
            for kk = 1:length(lift_bounds(:,1))
                all_trials{ii}(jj,:,kk) = FR{ii}(jj,...
                    lift_bounds(kk,1):lift_bounds(kk,2));
            end
        end
    end

    % Find the average FR for each unit across all trials.
    all_trials_averages = cell(num_probes,1);
    for ii = 1:num_probes
        all_trials_averages{ii} = mean(all_trials{ii},3);
    end

%% Plot CFA/RFA trial averages a series of figures
    
    % Determine how many figures are needed to graph the average FR for all
    % units in 4x4 subplots.
    num_plots = cell(num_probes,1);
    for ii = 1:num_probes
        num_plots{ii} = ceil(length(CFA_RFA_units{ii})/16);
    end

    % Plot every average FR in a series of figures
    CFA_RFA_FR = cell(num_probes,1);
    for ii = 1:num_probes
        CFA_RFA_FR{ii} = zeros(length(CFA_RFA_units{ii}),lift_duration);
    end
    time = lift_window(1):lift_window(2);
    for ii = 1:num_probes
        for jj = 1:num_plots{ii}
    %         fig = figure;
    %         fig.Units = 'inches';
    %         fig.Position = [3,2,10,8];
            for kk = 1:16
                % unit_ref is the index of unit and FR we are currently
                % iterating through. (This value may increase above
                % length(units) during the plotting, so a check is done below.
                unit_ref = (jj-1)*16+kk;
                if unit_ref <= length(CFA_RFA_units{ii})
                    %subplot(4,4,kk);
                    %text = sprintf('Probe %d: Unit %d',...
                        %ii,CFA_RFA_units{ii}(unit_ref));
                    % full_index is the index of the units array at the
                    % location of the current cortex unit iterate.
                    full_index = find(units{ii}==CFA_RFA_units{ii}(unit_ref));
                    CFA_RFA_FR{ii}(unit_ref,:) = ...
                        all_trials_averages{ii}(full_index,:);
    %                 % plot the individual CFA/RFA trial average on a subplot
    %                 plot(time,all_trials_averages{ii}(full_index,:));
    %                 title(text);
    %                 xline(0);
                end
                if kk == 16
                    %text = sprintf('Probe_%d_Part_%d.fig',ii,jj);
                    %savefig(text);
                end
            end
        end
    end

    
    
%% Take average of CFA RFA FR and plot
    
    if plot_FR_trialavg == 1
        figure(100);
        figure(101);
        for ii = 1:num_probes
            for jj = 1:length(CFA_RFA_FR{ii}(:,1))
                CFA_RFA_FR{ii}(jj,:) = CFA_RFA_FR{ii}(jj,:)-mean(CFA_RFA_FR{ii}(jj,1:50));
                if ii ==1
                    figure(100);
                    hold on
                    plot(CFA_RFA_FR{ii}(jj,:))
                    title('All CFA FR (Averaged Across Trials)')
                elseif ii == 2
                    figure(101);
                    hold on
                    plot(CFA_RFA_FR{ii}(jj,:))
                    title('All RFA FR (Averaged Across Trials)')
                end
            end
        end
    end
    avg_CFA_FR = mean(CFA_RFA_FR{1},1);
    avg_RFA_FR = mean(CFA_RFA_FR{2},1);
    avg_CFA_RFA_FR = [avg_CFA_FR; avg_RFA_FR];

    b = zeros(length(avg_CFA_FR),1,2);
    b(:,1,1) = .25*std(CFA_RFA_FR{1},0,1);
    b(:,1,2) = .25*std(CFA_RFA_FR{2},0,1);

    figure;
    boundedline(time,avg_CFA_RFA_FR,b);

    legend('CFA','RFA','Location','Northwest');
    title('Averaged CFA and RFA Neural Activity');
    xlabel('Reach Initiation at 0');
    xline(0);

%% Display EMG averages (muscles individually)

    % Time the EMG matrix so that only EMG in liftbounds windows is used. Will
    % create a (num_muslces x lift_duration x num_good_trials) matrix which can
    % then be trial averaged.

    all_EMG_trials = zeros(num_muscles,lift_duration,num_good_trials);
    for ii = 1:num_muscles
        for jj = 1:length(lift_bounds(:,1))
            all_EMG_trials(ii,:,jj) = EMG(ii,...
                lift_bounds(jj,1):lift_bounds(jj,2));
        end
    end

    % Find the average EMG for muscle across all trials.
    all_EMG_trials_averages = mean(all_EMG_trials,3);

    figure;
    for ii = 1:num_muscles
        if ii < 3
            plot(time,all_EMG_trials_averages)
            hold on
        else
            plot(time,all_EMG_trials_averages)
            hold on
        end
    end

    legend("Deltoid","Pectoralis","Bicep","Tricep","ECR","PL",'Location',...
        'Northwest')
    title('EMG trial averages')
    xlabel('Reach Initiation at 0');
    xline(0);

%% Plot Average CFA, RFA, and EMG

    figure;
    all_EMG_average = mean(all_EMG_trials_averages,1);
    yyaxis right
    plot(time,all_EMG_average,'k--');
    hold on
    yyaxis left
    plot(time,avg_CFA_FR,'b-');
    plot(time,avg_RFA_FR,'r-');
    legend('CFA','RFA','EMG','Location','Northwest');
    xline(0);
    title('Average CFA\RFA and EMG');
end