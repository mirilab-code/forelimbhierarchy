%% load data

% data = 'Z:\akiko\10052020_ss2_g0\preprocess\';
% addpath(data);
% addpath(genpath('Z:\Scripts\neuropixel'));
 %addpath('C:\Users\mirilab\Box\Miri Lab Shared\Scripts\Hannah\JCCG results\A thesis folder');
% addpath('Z:\Scripts\statistics');
% 
% data_path = data;
% home_path = 'C:\Users\mirilab\Box\Miri Lab Shared\Scripts\Mark\infotheory\';
% 
% files = dir([data_path '*.mat']);

% 
% for i=1:length(files)
%     load(files(i).name,'-mat');
% end
% cd(home_path);
%function [nullTE,allowed_shifts ]= TE_analysis (events, EMG)

addpath('C:\Users\mirilab\Box\Miri Lab Shared\Scripts\Hannah\JCCG results\A thesis folder');

fr_threshold = 1;   % 1 Hz??
nROIs = length(events);
duration = 0;
trains = {};
firingrates = {};
nUnits = zeros(1,nROIs);
[CFA, RFA]= events_to_cor(events, 1500);
trains{1}= CFA;
trains {2}= RFA;

for i=1:length(events)
    %trains{i} = events_to_train(events{i});
    %trains{i} = trains{i}(~cellfun(@isempty,trains{i}));
    fr = cellfun(@length,trains{i}) / max(cellfun(@max,trains{i})) * 1000; % the *1000 makes it in Hertz
    firingrates{i} = fr;
    
    % take out units with low firing rates
    bad_fr = find(fr < fr_threshold);
    trains{i}(bad_fr) = [];
    firingrates{i}(bad_fr) = [];
    
    nUnits(i) = length(trains{i});
    this_duration = max(cellfun(@max,trains{i}));
    if(this_duration > duration)
        duration = this_duration;
    end

end

disp('done!');
%%
TEMatrix = zeros(nUnits(1),nUnits(2));
binTrains = {};
delay = 30;
for t=1:nROIs
    binTrains{t} = [];
    for i=1:nUnits(t)
        disp([t i]);
        this_train = trains{t}{i};
        bt = train_to_binary(this_train,ceil(duration));
%         size(bt)
        binTrains{t} = [binTrains{t} ; bt];
    end
end

stacked = [binTrains{1} ; binTrains{2}];

totalN = sum(nUnits);
disp('done converting to binary trains!!');

%% get movement epochs
s = sum(EMG,1);
plot(s);
movement = movement_only(EMG,3007000,3012000);
% but we also want some period of time before movement when RFA and CFA are active
premove = 100;
starts = find(diff(movement)==1);
for i=1:length(starts)
    movement(starts(i)-premove:starts(i)) = 1;
end

% get just movement initiation -150 to +100
initiation = movement*0;
for i=1:length(starts)
    t = starts(i);
    initiation(t-150:t+100) = 1;
end

not_movement = ~movement;

%%
stacked_initiation = [];
slices = [];
for i=1:totalN
    disp([i totalN]);
    [new_binary,slices] = cattrain(initiation,stacked(i,:));
    stacked_initiation = [stacked_initiation; new_binary];
end
% so stacked initiation is the concatenated binary trains BUT NOT PADDED
% cattrains are padded in calculate_null_TE



%%  testing
%[observedTE,nullTE,allowed_shifts] = calculate_null_TE(stacked_initiation,slices,30,15,true);


%  testing 2
%[~,nullTE2,allowed_shifts2] = calculate_null_TE(stacked_initiation,slices,30,85,false);

%nullTE = cat(3,nullTE,nullTE2);
%allowed_shifts = cat(3,allowed_shifts,allowed_shifts2);
%end







%%