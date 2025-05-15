%% Data to Load: FULL EMG traces (from preprocessing), reach bounds
addpath('Z:\Scripts\boundedline');
addpath('Z:\Scripts\boundedline\Inpaint_nans');
addpath('Z:\Scripts\boundedline\boundedline');
%addpath('C:\Users\mirilab\OneDrive\Documents\MATLAB')
% load('success_failed_reaches.mat');

%Change str1 to appropriate path
%str1 = 'C:\Users\mirilab\OneDrive\Documents\MATLAB\EMG_reaching_analysis_Fig1/';
str1 = 'Z:\Reaching_EMG_analysis\EMG_dual_probe_Fig1\';
str2 = '\1000windowtrials';
str3 = 'Z:\Reaching_EMG_analysis\EMG_inactivation_Fig1\';


animal_list = {'a048';'a050';'ss2';'MA1';'a051';'MA2';'re7';'re9';'re8';'re12';'re10'};
%%
% Select wheher to sample reach data
% 1 = sample; 0 = no sample
sample = 1;

% for z = 1:length(animal_list) % Run all animals (remove comment on "end" at bottom of script)
%Select animal
z = 3;
animal = animal_list{z,1};

if strcmp(animal,'a048')
    lift_window = [-110,110];
    dates = ['03272020';'03302020';'03312020';'04012020';'04022020';'04032020'];

elseif strcmp(animal,'a050')
    lift_window = [-100,100];
    dates = ['09222020';'09232020'];

elseif strcmp(animal,'ss2')
    lift_window = [-130,130];
    dates = ['09302020';'10012020';'10022020';'10052020';'10062020'];

elseif strcmp(animal,'MA1')
    lift_window = [-150,150];
    dates = ['10242020';'10262020';'10272020';'10282020'];

elseif strcmp(animal,'a051')
    lift_window = [-300,300];
    dates = ['11102020';'11122020'];

elseif strcmp(animal,'MA2')
    lift_window = [-200,200];
    dates = ['11152020';'11172020'];

elseif strcmp(animal,'re7')
    lift_window = [-180,180];
    dates = ['10132021';'10142021';'10152021';'10182021';'10192021'];

elseif strcmp(animal,'re8')
    lift_window = [-100,100];
    dates = ['10272021';'10282021';'11012021'];
   
elseif strcmp(animal,'re9')
    lift_window = [-110,110];
    dates = ['11162021';'11172021';'11182021';'11192021'];
  
elseif strcmp(animal,'re12')
    lift_window = [-150,150];
    dates = ['12162021';'12172021'];

elseif strcmp(animal,'re10')
    lift_window = [-130,130];
    dates = ['01182022';'01192022';'01202022';'01212022'];
    
else
    disp('Error! Wrong date!')
end

dirs = [];
EMG_files = [];
reach_files = [];
spout_files = [];
for ii = 1:size(dates,1)
    date = dates(ii,:);
    if z<7
        dirs = [dirs;str1,date,'_reaching_analysis',str2];
    else
        dirs = [dirs;str3,date,'_reaching_analysis',str2];
    end
    EMG_files = [EMG_files;date,'_1000_EMG.mat'];
    reach_files = [reach_files;date,'_1000_reach_bounds_edit.mat'];
    spout_files = [spout_files;date,'_1000_spout_max_index_edit.mat'];
end
num_sessions = size(dates,1);
EMG_agg = cell(num_sessions,1);
reach_bounds_agg = cell(num_sessions,1);
spout_max_index_agg = cell(num_sessions,1);

for ii = 1:num_sessions
    cd(dirs(ii,:))
    load(EMG_files(ii,:));
    load(reach_files(ii,:));
    load(spout_files(ii,:));
    EMG_agg{ii} = EMG;
    reach_bounds_agg{ii} = reach_bounds_edit;
    spout_max_index_agg{ii} = spout_max_index_edit;
    clear EMG
    clear reach_bounds_edit
    clear spout_max_index_edit
end
mat_reach_bounds_agg = cell2mat(reach_bounds_agg);
reach_durations = mat_reach_bounds_agg(:,2) - mat_reach_bounds_agg(:,1);
figure;
hist(reach_durations);
disp(mean(reach_durations))
disp(median(reach_durations))

%%
animal_avg = cell(2,1);
b_all_trials = cell(2,1);
b = cell(2,1);
for align = 1:2
%align = 1; %% 1-aligned to movement 2-aligned to touch
% if align == 1
%     lift_window = [-150,150];  
% elseif align == 2
%     lift_window = [-150,150]; 
% end
  
full_time = cell(num_sessions,1);
num_good_trials = cell(num_sessions,1);
num_muscles = size(EMG_agg{1},1);
num_spouts = length(spout_max_index_agg{1});
EMGzscore = cell(num_sessions,1);
muscle_spout_EMG = cell(num_sessions,1);
b_all_trials{align} = cell(num_muscles,num_spouts);
for session = 1:num_sessions
    EMG = EMG_agg{session};
    reach_bounds_edit = reach_bounds_agg{session};
    full_time{session} = size(EMG,2);
    num_good_trials{session} = size(reach_bounds_edit,1);

    if (num_muscles == 4)
        muscles = {'Bicep','Tricep','ECR','PL'};
    elseif (num_muscles == 6)
        muscles = {'Deltoid','Pectoralis','Bicep','Tricep','ECR','PL'}; 
    end

    EMGzscore{session} = zeros(num_muscles,full_time{session});
    for ii = 1:num_muscles
        EMGzscore{session}(ii,:) = zscore(EMG(ii,:));
    end

    % We create windows -2000ms to +400ms surrounding the reach startclear
    lift_duration = lift_window(2)-lift_window(1)+1;
    lift_bounds = zeros(length(reach_bounds_edit(:,1)),2);

    % lift_bounds is a (number of successful trials for all spouts x 2) matrix
    % that contains all of the lift windows for all accepted trials.
    for ii = 1:length(reach_bounds_edit(:,1))
        lift_bounds(ii,:) = [reach_bounds_edit(ii,align)+lift_window(1),...
           reach_bounds_edit(ii,align)+lift_window(2)];
    end

    all_EMG_trials_z = zeros(num_muscles,lift_duration,num_good_trials{session});
    for ii = 1:num_muscles
        for jj = 1:length(lift_bounds(:,1))
            all_EMG_trials_z(ii,:,jj) = EMGzscore{session}(ii,...
                lift_bounds(jj,1):lift_bounds(jj,2));
        end
    end
    
    spout_max_index = spout_max_index_agg{session};
    time = lift_window(1):lift_window(2);
 
    EMG_spout_avg = zeros(num_muscles,lift_duration,num_spouts);
    running = 1;
    all_avg = mean(all_EMG_trials_z,3);

    %% plot on different spout
        
    figure;
    for ii = 1:num_spouts
        if spout_max_index{ii} ~= 0
            EMG_spout_avg(:,:,ii) = ...
                mean(all_EMG_trials_z(:,:,running:spout_max_index{ii}),3);
            %running = spout_max_index{ii}+1;
            for jj = 1:num_muscles
%                 b_all_trials{align}{jj,ii} = ...
%                     [b_all_trials{align}{jj,ii};all_EMG_trials_z(:,:,running:spout_max_index{ii})];
            a = all_EMG_trials_z(jj,:,running:spout_max_index{ii});
            b_all_trials{align}{jj,ii} = ...
                     cat(3,b_all_trials{align}{jj,ii},a);
            end
            running = spout_max_index{ii}+1;
        end
        muscle_spout_EMG{session} = EMG_spout_avg;
        subplot(num_muscles/2,2,ii);
        hold on
        for jj = 1:num_muscles
            plot(time,EMG_spout_avg(jj,:,ii));
        end
        str = sprintf('Spout %d',ii);
        title(str)
        if (num_muscles == 4)
            legend('Bicep','Tricep','ECR','PL','Location','Southeast');
        elseif (num_muscles == 6)
            legend('Deltoid','Pectoralis','Bicep','Tricep','ECR','PL',...
            'Location','Southeast'); 
        end  
        xlabel('Time (ms)')
        ylabel('EMG (zscored)')
    end
    str = sprintf('%s: Trial Average Muscle EMG',dates(session,:));
    sgtitle(str);
    
    figure;
    for ii = 1:num_muscles
        for jj = 1:num_spouts
            subplot(num_muscles/2,2,ii);
            hold on
            plot(time,EMG_spout_avg(ii,:,jj));
        end
        title(muscles{ii})
        legend('spout 1','spout 2','spout 3','spout 4','Location','Southeast');
        xlabel('Time (ms)')
        ylabel('EMG (zscored)')
    end
    str = sprintf('%s: Trial Average Muscle EMG',dates(session,:));
    sgtitle(str);
    
    figure;
    hold on
    for ii = 1:num_spouts
        plot(time,mean(EMG_spout_avg(:,:,ii),1));
    end
    legend('spout 1','spout 2','spout 3','spout 4','Location','Southeast');
    str = sprintf('%s: Trial and Muscle Averaged EMG',dates(session,:));
    title(str)
    xlabel('Time (ms)')
    ylabel('EMG (zscored)')

    figure;
    hold on
    for ii = 1:num_muscles
        plot(time,mean(EMG_spout_avg(ii,:,:),3));
    end
    if (num_muscles == 4)
       legend('Bicep','Tricep','ECR','PL','Location','Southeast');
    elseif (num_muscles == 6)
       legend('Deltoid','Pectoralis','Bicep','Tricep','ECR','PL',...
           'Location','Southeast'); 
    end
    str = sprintf('%s: Trial and Spout Averaged EMG',dates(session,:));
    title(str)
    xlabel('Time (ms)')
    ylabel('EMG (zscored)')
end


total_num_trials = sum(cell2mat(num_good_trials));
animal_avg{align} = zeros(num_muscles,lift_duration,num_spouts);
b{align} = zeros(num_muscles,lift_duration,num_spouts);

for ii = 1:num_muscles
    for jj = 1:num_spouts
        running_traces = [];
        running = 0;
        for kk = 1:num_sessions
            running = running + (num_good_trials{kk}/total_num_trials)*...
                muscle_spout_EMG{kk}(ii,:,jj);
        end
        animal_avg{align}(ii,:,jj) = running;
    end
end
end

%%

ylims = zeros(num_muscles,2);
for ii = 1:num_muscles
    temp = cell2mat(animal_avg');
    ylims(ii,2) = max(temp(ii,:))+0.1*max(temp(ii,:));
    ylims(ii,1) = min(temp(ii,:))-0.1*max(temp(ii,:));
end

%%
%Determine number of trials for each spout
spout_num_trials_total = [size(b_all_trials{1,1}{1,1},3), size(b_all_trials{1,1}{1,2},3), ...
    size(b_all_trials{1,1}{1,3},3), size(b_all_trials{1,1}{1,4},3)];


% Concatenate b_all_trials for plotting
b_all_trials_plot_2 = cell(2,1);
for iii = 1:2 % align num
    for ii = 1:4 % spout num
        b_all_trials_temp = [];
        for i = 1:num_muscles % muscle num
            b_all_trials_temp = [b_all_trials_temp; mean(b_all_trials{iii,1}{i,ii},3)];
        end

        b_all_trials_plot_2{iii,1}(:,:,ii) = b_all_trials_temp;
    end
end
%%
close all
figure;
box off
hold on
for ii = 1:2*num_muscles*num_spouts
    subplot(num_muscles,2*num_spouts,ii)
    if mod(ii,2) == 1
        ii = (ii+1)/2;
        align = 1;
        color = 'k';
    else
        ii = ii/2;
        align = 2;
        color = 'g';
    end
    q = floor(ii/num_spouts-0.001);
    r = mod(ii,num_spouts);
    if r == 0
        r = 4;
    end
    trace = b_all_trials_plot_2{align}(q+1,:,r);
    time = -1*((size(b_all_trials_plot_2{align},2)-1)/2):(size(b_all_trials_plot_2{align},2)-1)/2;
    b = std(b_all_trials{align}{q+1,r},0,3)/sqrt(spout_num_trials_total(r));
    %plot(time,trace,color);
    boundedline(time,trace,b,color)
    xlim(lift_window)
    ylim(ylims(q+1,:))
    if r == 1 && align == 1
        ylabel(muscles{q+1})
    end
    if q+1 == num_muscles && align == 1
        str = sprintf('Spout %d',r);
        xlabel(str)
    end
    
    if ii == 1
        title(mat2str(spout_num_trials_total(1,1)));
    elseif ii == 2
        title(mat2str(spout_num_trials_total(1,2)));
    elseif ii == 3
        title(mat2str(spout_num_trials_total(1,3)));
    elseif ii == 4
        title(mat2str(spout_num_trials_total(1,4)));
    end
end
str = sprintf('%s: black - movement onset, green - beam break',animal);
sgtitle(str)
cd 'C:\Users\mirilab\OneDrive\Documents\MATLAB';
saveas(gcf,['muscle_avg_',animal_list{z,1},'.png'])
%% Sample Data
if sample == 1
% Normalize reaches to every spout to the lowest reached spout
[small_reach, small_reach_index] = min(spout_num_trials_total);
prune_index_s1 = randperm(size(b_all_trials{1,1}{1,1},3), small_reach);
prune_index_s2 = randperm(size(b_all_trials{1,1}{1,2},3), small_reach);
prune_index_s3 = randperm(size(b_all_trials{1,1}{1,3},3), small_reach);
prune_index_s4 = randperm(size(b_all_trials{1,1}{1,4},3), small_reach);
prune_index = [prune_index_s1; prune_index_s2; prune_index_s3; prune_index_s4]; % rows are spout number

b_all_trials_new = cell(2,1);
for iii = 1:2 % align num
    for ii = 1:4 % spout num
        for i = 1:num_muscles % muscle num
            b_all_trials_temp = [];
            for iiii = prune_index(ii,:) % reach num
                b_all_trials_temp = cat(3, b_all_trials_temp, b_all_trials{iii,1}{i,ii}(:,:,iiii));
            end
            
            b_all_trials_new{iii,1}{i,ii} = b_all_trials_temp;
        end
    end
end



%%
% Concatenate b_all_trials for plotting
b_all_trials_plot = cell(2,1);
for iii = 1:2 % align num
    for ii = 1:4 % spout num
        b_all_trials_temp = [];
        for i = 1:num_muscles % muscle num
            b_all_trials_temp = [b_all_trials_temp; mean(b_all_trials_new{iii,1}{i,ii},3)];
        end

        b_all_trials_plot{iii,1}(:,:,ii) = b_all_trials_temp;
    end
end

%%
figure;
box off
hold on
for ii = 1:2*num_muscles*num_spouts
    subplot(num_muscles,2*num_spouts,ii)
    if mod(ii,2) == 1
        ii = (ii+1)/2;
        align = 1;
        color = 'k';
    else
        ii = ii/2;
        align = 2;
        color = 'g';
    end
    q = floor(ii/num_spouts-0.001);
    r = mod(ii,num_spouts);
    if r == 0
        r = 4;
    end
    small_reach_sem = [small_reach, small_reach, small_reach, small_reach];
    trace = b_all_trials_plot{align}(q+1,:,r);
    time = -1*((size(b_all_trials_plot{align},2)-1)/2):(size(b_all_trials_plot{align},2)-1)/2;
    b = std(b_all_trials_new{align}{q+1,r},0,3)/sqrt(small_reach_sem(r));
    %plot(time,trace,color);
    boundedline(time,trace,b,color)
    xlim(lift_window)
    ylim([-0.7, 2.6])
%     ylim(ylims(q+1,:))
    if r == 1 && align == 1
        ylabel(muscles{q+1})
    end
    if q+1 == num_muscles && align == 1
        str = sprintf('Spout %d',r);
        xlabel(str)
    end

%     if ii == 1 || ii == 2 || ii == 3 || ii == 4
%         title(mat2str(small_reach));
%     end

end
str = sprintf('%s: black - movement onset, green - beam break',[animal,' sampled']);
sgtitle(str)
cd 'C:\Users\mirilab\OneDrive\Documents\MATLAB';
saveas(gcf,['muscle_avg_samp_',animal_list{z,1},'.png'])
end
% end

