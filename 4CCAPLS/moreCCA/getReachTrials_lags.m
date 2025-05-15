%% CCA 


disp('loading things...');

addpath('C:\Users\mirilab\Documents\Mark\util');

folders = dir('Z:\akiko');
datafolders = struct2cell(folders);
foldernames = datafolders(1,3:end-10)';

opts = detectImportOptions('trialstats.csv');
opts = setvartype(opts, 'Date', 'char');
stats = readtable('trialstats.csv', opts);
dates = stats.Date(1:end-2);
animals = stats.Animal(1:end-2);
%%
animal = 'MA2';
% lags = -30:30;    % so lag is how much CFA is moved, see below
lags = 5:30;
aidx = find(contains(animals,animal));

fidx = find(contains(foldernames,animal));
exps = foldernames(fidx);

tic
for lag=lags
    CFA = [];
    RFA = [];
    disp(lag);
    totalCells1 = 0;
    totalCells2 = 0;
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


        disp('converting to firing rates...');
        nProbes = length(events);
        duration = size(EMG,2);
        trains = {};
        FR = {};
        for i=1:nProbes
            trains{i} = events_to_train(events{i});
            FR{i} = train_to_FR(trains{i});
            FR{i} = FR{i}(:,1:duration);
        end
        nCells1 = length(trains{1});
        nCells2 = length(trains{2});
        totalCells1 = totalCells1 + nCells1;
        totalCells2 = totalCells2 + nCells2;
        fprintf('%d in CFA and %d in RFA \n', nCells1, nCells2);
        fprintf('for a total of %d in CFA and %d in RFA \n', totalCells1, totalCells2);



        window = -99:100;

        oneexpCFA = {};
        oneexpRFA = {};
        disp('separating reaches and grasps...');
        for s=1:length(spout_bounds)
        %     disp(s);
            trials = spout_bounds{s};
            nTrials = size(trials,1);
            oneexpCFA{s} = zeros(nCells1,2*length(window),nTrials);
            oneexpRFA{s} = zeros(nCells2,2*length(window),nTrials);
            for t=1:nTrials
        %         disp([s t]);
                reach = trials(t,1) + window;
                grasp = trials(t,2) + window;
                traceCFA = [FR{1}(:,reach+lag) FR{1}(:,grasp+lag)];     % so lag is how much CFA is moved
                traceRFA = [FR{2}(:,reach) FR{2}(:,grasp)];
                oneexpCFA{s}(:,:,t) = traceCFA;
                oneexpRFA{s}(:,:,t) = traceRFA;
            end
            oneexpCFA{s} = mean(oneexpCFA{s},3);
            oneexpRFA{s} = mean(oneexpRFA{s},3);
        end

        stackedCFA = {oneexpCFA{1} oneexpCFA{2} oneexpCFA{3} oneexpCFA{4}};
        stackedCFA = cat(2,stackedCFA{:});
        stackedRFA = {oneexpRFA{1} oneexpRFA{2} oneexpRFA{3} oneexpRFA{4}};
        stackedRFA = cat(2,stackedRFA{:});

        disp('concatenating experiments...');
        CFA = [CFA; stackedCFA];
        RFA = [RFA; stackedRFA];
    end
    disp('done!');
    %%
    disp('saving...');
    eval([animal ' = struct']);
    eval([animal '.CFA = CFA']);
    eval([animal '.RFA = RFA']);

    savename = sprintf('output/%s_%d.mat',animal,lag);
    save(savename,animal);
    
end
disp('done!!');
toc
%%







%%