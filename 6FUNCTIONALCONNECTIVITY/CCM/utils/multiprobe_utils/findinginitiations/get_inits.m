function [init_times] = diya_get_inits

clear all
close all
addpath(genpath('Z:/code'));
%%

disp('go to preprocess folder')
df = uigetdir('Z:\multiprobe_climbing_data');
load_all_from_dir(df)

%%
encoder_raw = analogin(2,:);
plot(encoder_raw)

%%
immobility_min_duration = 0.4; %minimal duration of "immobility" periods (sec)
immobility_min_duration_before_climbing = 1; %minimal duration of the immobility period preceding "climbing" bouts (sec)
min_climbing_angle = 10; %minimum angular distance travelled during a climbing bout (°)
[climbing_trials,encoder_filtered,isclimbing]= ...
        pre_processing_wheel(encoder_raw,immobility_min_duration,immobility_min_duration_before_climbing,min_climbing_angle);

init_times = [climbing_trials.onset];


%% Save
cd('Z:\Sarah\multiprobe_DLAG\')

save_path = uigetdir('', 'Save Location');
cd(save_path)

save('climbing_trials.mat', 'climbing_trials');






























%%