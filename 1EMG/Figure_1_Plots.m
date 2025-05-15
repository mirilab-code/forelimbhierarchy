%% Figure 1 Plots

%% Time Course Schematic
addpath('C:\Users\mirilab\OneDrive\Documents\MATLAB')
load('example_EMG');

emg = zeros(2500,1);
emg_y = [example_EMG(3,1:411), example_EMG(3,1:400), example_EMG(3,1:300), example_EMG(3,1:1189), flip(example_EMG(3,990:1189))]';
emg(1:length(emg_y),1) = emg_y;

led = zeros(2500,1);
led_y = (0.5*ones(2100,1));
led(300:300+length(led_y)-1,1) = led_y;

tone = zeros(2500,1);
tone_y = (0.5*ones(100,1));
tone(1500:1500+length(tone_y)-1,1) = tone_y;

solenoid = zeros(2500,1);
solenoid_y = (0.5*ones(100,1));
solenoid(1600:1600+length(solenoid_y)-1,1) = solenoid_y;

beam = zeros(2500,1);
beam_y = 0.5*ones(1,1);
beam(1700:1700+length(beam_y)-1,1) = beam_y;

figure;
plot(led(1:2500,1)+5,'color','k','LineWidth',1);
hold on
plot(tone(1:2500,1)+4,'color','k','LineWidth',1);
hold on
plot(solenoid(1:2500,1)+3,'color','k','LineWidth',1);
hold on
plot(beam(1:2500,1)+2,'color','k','LineWidth',1);
hold on
plot((emg(1:2500,1)/100)+1,'color','k','LineWidth',1);
set(gca,'ytick',[0:10],'yticklabel',{'';'EMG';'Beam Break';'Solenoid';'Reward Tone';'LED';'';'';'';'';''});
% set(gca,'xtick',[]);
% set(gca,'tickdir','none');
xlabel('Time (ms)');

%% Trial Averaged Muscles for All Spouts

%% Behavior Performance

% Add paths to needed functions in this script
% (These paths need to be edited to the proper location of these functions)
opengl software
home_path = 'Z:\Scripts';
addpath(home_path);
addpath('Z:\neuropixel');
neur_paths = genpath('Z:\Scripts');
addpath(neur_paths);
addpath('C:\Users\mirilab\OneDrive\Documents\MATLAB\util');

% Global Variables
% Recording Identifiers
dates = ['03272020';'03302020';'03312020';'04012020';'04022020';...
    '04032020';'09222020';'09232020';'09302020';'10012020';'10022020';...
    '10052020';'10062020';'10242020';'10262020';'10272020';'10282020';...
    '11102020';'11122020';'11152020';'11172020'];
% dates = ['03302020'];

animals = {'a048';'a048';'a048';'a048';'a048';'a048';'a050';'a050';...
    'ss2';'ss2';'ss2';'ss2';'ss2';'MA1';'MA1';'MA1';'MA1';...
    'a051';'a051';'MA2';'MA2'};
% animals = {'a048'};

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



% Create List Of # of Successful and Failed Reaches (Excluding Non-Attempts)
allowable_error_range = 50;

error_best_spout = cell(length(allowable_error_range)+1,length(unique(animals)));
error_second_spout = cell(length(allowable_error_range)+1,length(unique(animals)));
error_best_spout(1,:) = unique(animals,'stable');
error_second_spout(1,:) = unique(animals,'stable');

counter_error = 2;
for error_new = allowable_error_range

success_failed_reaches = cell(4,length(dates));
for dates_animals = 1:size(dates,1)
    date = dates(dates_animals,:);
    animal = animals{dates_animals};
    fprintf('recording %s, animal %s\n',date,animal);

    success_failed_reaches{1,dates_animals} = date;
    success_failed_reaches{2,dates_animals} = animal;

    % Load Data
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

    % Generate reward and units cells

    % Create reward cell containing the times when the reward sound is
    % detected (this is for both successful and unsuccesful trials).
    rw = get_reward_times(digital, [reward1 reward2 reward3 reward4]);

    % correct for the reward sound by subtracting 100
    rw1 = rw{reward1}-100;
    rw2 = rw{reward2}-100;
    rw3 = rw{reward3}-100;
    rw4 = rw{reward4}-100;
    rw_corr = {rw1,rw2,rw3,rw4};

    % Getting reach data
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
    allowable_error = error_new;
    success_bounds = find_success(num_spouts,rw_corr,touch_channels,...
        digital(suction,:),allowable_error);
    
    % find_failed() determines the number of failed reaches to each spout
    failed_reaches = find_failed(num_spouts,rw_corr,touch_channels,...
        digital(suction,:),allowable_error);

    % Populate success and failed reaches per animal per spout
    success_failed_reaches{3,dates_animals} = [length(success_bounds{1,1}), length(success_bounds{1,2}),...
        length(success_bounds{1,3}), length(success_bounds{1,4})];
    success_failed_reaches{4,dates_animals} = [failed_reaches{1,1}, failed_reaches{1,2},...
        failed_reaches{1,3},failed_reaches{1,4}];
end

% Concatenate all success and all failures between mice per spout
success_total = [];
failed_total = [];
for i = 1:length(success_failed_reaches)
    success_total(i,1:4) = success_failed_reaches{3,i};
    failed_total(i,1:4) = success_failed_reaches{4,i};
end

success_percent = (success_total./(success_total + failed_total)) * 100;
time = [1:size(success_percent,1)];
% plot individual spout and recording session success rate
figure();
subplot(4,1,1)
bar(success_percent(:,1));
title('Success Rate: Spout 1')
xlabel('Recording Session')
ylabel('Success Rate (%)')
hold on

subplot(4,1,2)
bar(success_percent(:,2));
title('Success Rate: Spout 2')
xlabel('Recording Session')
ylabel('Success Rate (%)')
hold on

subplot(4,1,3)
bar(success_percent(:,3));
title('Success Rate: Spout 3')
xlabel('Recording Session')
ylabel('Success Rate (%)')
hold on

subplot (4,1,4)
bar(success_percent(:,4));
title('Success Rate: Spout 4')
xlabel('Recording Session')
ylabel('Success Rate (%)')
hold off

% plot averaged recording session per spout success rate
time = [1:size(success_percent,2)];
for i = 1:size(success_percent,2)
    success_percent_avg_spout(1,i) = mean(success_percent(:,i));
end

figure();
bar(success_percent_avg_spout);
title('Avg Success Rate Per Spout')
xlabel('Spout')
ylabel('Avg Success Rate (%)')

% plot averaged spout success per session
time = [1:size(success_percent,1)];
for i = 1:size(success_percent,1)
    success_percent_avg_session(i) = mean(success_percent(i,:));
end

figure();
bar(success_percent_avg_session);
title('Avg Success Rate Per Recoding Session')
xlabel('Recording Session')
ylabel('Avg Success Rate (%)')

% weighted average per animal all spouts    
[c,ia,ic] = unique(success_failed_reaches(2,:),'stable');
success_rate_weighted_avg = cell(2,length(ia));
for i = 1:length(ia)
    success_rate_weighted_avg{1,i} = success_failed_reaches {2,ia(i)};
    if i == length(ia)
        success_rate_weighted_avg{2,i} = (sum(success_total(ia(i):length(success_failed_reaches),:),'all')/...
            (sum(success_total(ia(i):length(success_failed_reaches),:),'all')...
            + sum(failed_total(ia(i):length(success_failed_reaches),:),'all'))) * 100;
    else
        success_rate_weighted_avg{2,i} = (sum(success_total(ia(i):ia(i+1)-1,:),'all')/...
            (sum(success_total(ia(i):ia(i+1)-1,:),'all') + sum(failed_total(ia(i):ia(i+1)-1,:),'all'))) * 100;
    end
end

% average best and second best spout per animal
success_rate_best_spout = cell(2,length(ia));
success_rate_second_spout = cell(2,length(ia));
for i = 1:length(ia)
    success_rate_best_spout{1,i} = success_failed_reaches{2,ia(i)};
    success_rate_second_spout{1,i} = success_failed_reaches{2,ia(i)};
    success_rate_animal_spout = [];
    if i == length(ia)
        success_rate_animal_spout = (sum(success_total(ia(i):length(success_failed_reaches),:))./...
            (sum(success_total(ia(i):length(success_failed_reaches),:))...
            + sum(failed_total(ia(i):length(success_failed_reaches),:)))) .* 100;
    else
        success_rate_animal_spout = (sum(success_total(ia(i):ia(i+1)-1,:))./...
            (sum(success_total(ia(i):ia(i+1)-1,:)) + sum(failed_total(ia(i):ia(i+1)-1,:)))) .* 100;
    end
    success_rate_best_spout{2,i} = max(success_rate_animal_spout);
    second = sort(success_rate_animal_spout, 'descend');
    success_rate_second_spout{2,i} = second(2);
end

    % plot
    xscale_best_spout =  ones(1,length(success_rate_best_spout));
    xscale_weighted_avg = 3 .* ones(1,length(success_rate_weighted_avg));
    xscale_second_spout = 2 .* ones(1,length(success_rate_weighted_avg));
    figure();
    % plot scatters
    scatter(xscale_best_spout,cell2mat(success_rate_best_spout(2,:)),75,'filled');
    hold on
    scatter(xscale_weighted_avg,cell2mat(success_rate_weighted_avg(2,:)),75,'filled');
    hold on
    scatter(xscale_second_spout,cell2mat(success_rate_second_spout(2,:)),75,'filled');
    % plot means
    plot(1,mean(cell2mat(success_rate_best_spout(2,:))),'k_','MarkerSize',12,'LineWidth',2);
    hold on
    plot(3,mean(cell2mat(success_rate_weighted_avg(2,:))),'k_','MarkerSize',12,'LineWidth',2);
    hold on
    plot(2,mean(cell2mat(success_rate_second_spout(2,:))),'k_','MarkerSize',12,'LineWidth',2);
    % plot connected lines
    for i = 1:length(success_rate_best_spout(2,:))
        plot([1, 2, 3], [cell2mat(success_rate_best_spout(2,i)), cell2mat(success_rate_second_spout(2,i)),...
            cell2mat(success_rate_weighted_avg(2,i))],'k','LineWidth',1);
        hold on
    end
    x_label = {'Best Spout','Second Best Spout','All Spouts'};
    set(gca,'XTick',1:3,'XTickLabel',x_label, 'Xlimitmethod', 'padded', 'Tickdir', 'out');
    ylim([0, 100]);
    title('Avg Success Rate Per Animal');
    ylabel('Success Rate (%)');

% Create vector of success rate of all mice on best and second best spout
% per session

[c,ia,ic] = unique(success_failed_reaches(2,:),'stable');
success_rate_best_spout = cell(2,length(ia));
success_rate_second_spout = cell(2,length(ia));
for i = 1:length(ia)
    success_rate_best_spout{1,i} = success_failed_reaches{2,ia(i)};
    success_rate_second_spout{1,i} = success_failed_reaches{2,ia(i)};
    success_rate_animal_spout = [];
    success_rate_best_spout{2,i} = [];
    success_rate_second_spout{2,i} = [];
    if i == length(ia)
        counter = 1;
        for ii = ia(i):length(success_failed_reaches)
            success_rate_animal_spout = success_total(ii,:) ./...
                (success_total(ii,:) + failed_total(ii,:)) .* 100;

            [B, I] = maxk(success_rate_animal_spout,2);
            success_rate_best_spout{2,i}(counter,1) = B(1,1);
            success_rate_best_spout{3,i}(counter,1) = I(1,1);

            success_rate_second_spout{2,i}(counter,1) = B(1,2);
            success_rate_second_spout{3,i}(counter,1) = I(1,2);


            counter = counter + 1;
        end
    else
        counter = 1;
        for ii = ia(i):ia(i+1)-1
            success_rate_animal_spout = success_total(ii,:)./...
                (success_total(ii,:) + failed_total(ii,:)) .* 100;

            [B, I] = maxk(success_rate_animal_spout,2);
            success_rate_best_spout{2,i}(counter,1) = B(1,1);
            success_rate_best_spout{3,i}(counter,1) = I(1,1);

            success_rate_second_spout{2,i}(counter,1) = B(1,2);
            success_rate_second_spout{3,i}(counter,1) = I(1,2);

            counter = counter + 1;
        end
    end
end

% Plot best spout and second best spout per session per mouse
figure();
for i = 1:length(ia)
    subplot(length(ia),1,i);
    b = bar([success_rate_best_spout{2,i}, success_rate_second_spout{2,i}]);
    title(['Success Rate ', success_rate_best_spout{1,i}]);
    ylabel('Success Rate (%)');
    xlabel('Recording Session');
    ylim([1, 125]);
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(round(b(2).YData));
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    hold on
end
legend('Best Spout','Second Best Spout', 'Location','eastoutside');

error_best_spout(counter_error,:) = success_rate_best_spout(2,:);
error_second_spout(counter_error,:) = success_rate_second_spout(2,:);
counter_error = counter_error + 1;

disp(['Finished with Error: ', mat2str(error_new)]);
end

% Plot correlation between success rate and allowable error for best spout

% Organize success rates by session
error_session_best = [];
error_session_line_best = [];
for mice = 1:size(error_best_spout,2)
    counter_session = 0;
    session_line = [];
    for session = 2:size(error_best_spout,1)
        error_session_best = [error_session_best,[]; error_best_spout{session,mice}, counter_session .* ones(length(error_best_spout{session,mice}),1)];
        counter_session = counter_session + 10;

        session_line = [session_line, error_best_spout{session,mice}];
    end
    error_session_line_best = [error_session_line_best; session_line];
end

xaxis = allowable_error_range;
cm = jet(size(error_session_line_best,1));
figure();
% plot scatter of all success rates per error
scatter(error_session_best(:,2),error_session_best(:,1));
hold on
% h1 = lsline;
% h1.Color = 'k';
% h1.LineWidth = 1;
% hold on
yaxis = mean(error_session_line_best,1);
plot(xaxis, yaxis, 'color', 'k','LineWidth', 1);
hold on
% plot lines connecting same sessions
for i = 1:size(error_session_line_best,1)
    plot(xaxis, error_session_line_best(i,:),'color',cm(i,:),'marker','o','MarkerFaceColor',cm(i,:));
    hold on
end
title('Correlation Between Success Rate and Error for Best Spout');
ylabel('Success Rate (%)');
xlabel('Allowable Error');

% Plot correlation between success rate and allowable error for second best spout

% Organize success rates by session
error_session_second = [];
error_session_line_second = [];
for mice = 1:size(error_second_spout,2)
    counter_session = 0;
    session_line = [];
    for session = 2:size(error_second_spout,1)
        error_session_second = [error_session_second,[]; error_second_spout{session,mice}, counter_session .* ones(length(error_second_spout{session,mice}),1)];
        counter_session = counter_session + 10;

        session_line = [session_line, error_second_spout{session,mice}];
    end
    error_session_line_second = [error_session_line_second; session_line];
end

xaxis = allowable_error_range;
cm = jet(size(error_session_line_second,1));
figure();
% plot scatter of all success rates per error
scatter(error_session_second(:,2),error_session_second(:,1));
hold on
% h1 = lsline;
% h1.Color = 'k';
% h1.LineWidth = 1;
% hold on
yaxis = mean(error_session_line_second,1);
plot(xaxis, yaxis, 'color', 'k','LineWidth', 1);
hold on
% plot lines connecting same sessions
for i = 1:size(error_session_line_second,1)
    plot(xaxis, error_session_line_second(i,:),'color',cm(i,:),'marker','o','MarkerFaceColor',cm(i,:));
    hold on
end
title('Correlation Between Success Rate and Error for Second Best Spout');
ylabel('Success Rate (%)');
xlabel('Allowable Error');

concatenate_best = [];
concatenate_second = [];
for i = 1:6
    concatenate_best = [concatenate_best; success_rate_best_spout{3,i}];
    concatenate_second = [concatenate_second; success_rate_second_spout{3,i}]
end

figure();
histogram(concatenate_best);

figure();
histogram(concatenate_second);



























