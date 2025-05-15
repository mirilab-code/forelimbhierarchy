%% make spout trial plots

%% CCA 
CFA = [];
RFA = [];
full_EMG = [];

disp('loading things...');

% addpath(genpath('Z:\Scripts'));

folders = dir('Z:\akiko');
datafolders = struct2cell(folders);
foldernames = datafolders(1,3:end-10)';

opts = detectImportOptions('trialstats.csv');
opts = setvartype(opts, 'Date', 'char');
stats = readtable('trialstats.csv', opts);
dates = stats.Date(1:end-2);
animals = stats.Animal(1:end-2);
%%
animal = 'ss2';
aidx = find(contains(animals,animal));

fidx = find(contains(foldernames,animal));
exps = foldernames(fidx);

totalCells1 = 0;
totalCells2 = 0;
numChannels = 4;

for d=1:length(exps)
    disp(exps{d});
    disp('loading things...');
    datapath = fullfile('Z:\akiko\',exps{d},'\preprocess_with_acg');
    files = fullfile(datapath,'*.mat');
    matFiles = dir(files);
    for i = 1:length(matFiles)
        baseFileName = fullfile(datapath, matFiles(i).name);
        load(baseFileName);
    end
    
    disp('getting trials bounds...');
    statsIdx = aidx(d);
    spoutNs = [stats.Spout1Used(statsIdx) stats.Spout2Used(statsIdx) stats.Spout3Used(statsIdx) stats.Spout4Used(statsIdx)];
    disp(spoutNs);
    spoutIdxs = cumsum(spoutNs);

    spout_bounds = {};
    spout_bounds{1} = reach_bounds_edit(1:spoutIdxs(1),:);
    spout_bounds{2} = reach_bounds_edit(spoutIdxs(1)+1:spoutIdxs(2),:);
    spout_bounds{3} = reach_bounds_edit(spoutIdxs(2)+1:spoutIdxs(3),:);
    spout_bounds{4} = reach_bounds_edit(spoutIdxs(3)+1:spoutIdxs(4),:);
    window = -199:200;

    oneexpCFA = {};
    oneexpRFA = {};
    oneexpEMG = {};
    disp('separating reaches and grasps...');
    for s=1:length(spout_bounds)
        trials = spout_bounds{s};
        nTrials = size(trials,1);
        oneexpEMG{s} = zeros(numChannels,2*length(window),nTrials);
        for t=1:nTrials
            reach = trials(t,1) + window;
            grasp = trials(t,2) + window;
            traceEMG = [EMG(:,reach) EMG(:,grasp)];
            oneexpEMG{s}(:,:,t) = traceEMG;
        end
        oneexpEMG{s} = mean(oneexpEMG{s},3);
    end
    stackedEMG = {oneexpEMG{1} oneexpEMG{2} oneexpEMG{3} oneexpEMG{4}};
    stackedEMG = cat(2,stackedEMG{:});
    disp('concatenating experiments...');
   
    full_EMG = cat(3,full_EMG,stackedEMG);
end
EMGtrials = mean(full_EMG,3);
disp('done!');


%%
[~,EMGscore,EMGlatent] = pca(EMGtrials');
plot(EMGscore(:,1));

%%
figure;
for i=1:8
    subplot(1,8,i)
    x = 400 * (i-1) +1;
    trc = EMGscore(x:x+399,1);
    plot(trc)
end
















%%