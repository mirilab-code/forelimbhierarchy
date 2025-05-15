%% make spout trial plots

%% CCA 
CFA = [];
RFA = [];
full_EMG = [];

disp('loading things...');

addpath('C:\Users\mirilab\Documents\Mark\util');
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
animal = 'a048';
aidx = find(contains(animals,animal));

fidx = find(contains(foldernames,animal));
exps = foldernames(fidx);

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
    nMuscles = size(EMG,1);
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
    oneexpEMG = {};
    disp('separating reaches and grasps...');
    for s=1:length(spout_bounds)
    %     disp(s);
        trials = spout_bounds{s};
        nTrials = size(trials,1);
        oneexpCFA{s} = zeros(nCells1,2*length(window),nTrials);
        oneexpRFA{s} = zeros(nCells2,2*length(window),nTrials);
        oneexpEMG{s} = zeros(nMuscles,2*length(window),nTrials);
        for t=1:nTrials
    %         disp([s t]);
            reach = trials(t,1) + window;
            grasp = trials(t,2) + window;
            traceCFA = [FR{1}(:,reach) FR{1}(:,grasp)];    
            traceRFA = [FR{2}(:,reach) FR{2}(:,grasp)];
            oneexpCFA{s}(:,:,t) = traceCFA;
            oneexpRFA{s}(:,:,t) = traceRFA;
            
            traceEMG = [EMG(:,reach) EMG(:,grasp)];
            oneexpEMG{s}(:,:,t) = traceEMG;
        end
        oneexpCFA{s} = mean(oneexpCFA{s},3);
        oneexpRFA{s} = mean(oneexpRFA{s},3);
        oneexpEMG{s} = mean(oneexpEMG{s},3);

    end

    stackedCFA = {oneexpCFA{1} oneexpCFA{2} oneexpCFA{3} oneexpCFA{4}};
    stackedCFA = cat(2,stackedCFA{:});
    stackedRFA = {oneexpRFA{1} oneexpRFA{2} oneexpRFA{3} oneexpRFA{4}};
    stackedRFA = cat(2,stackedRFA{:});
    stackedEMG = {oneexpEMG{1} oneexpEMG{2} oneexpEMG{3} oneexpEMG{4}};
    stackedEMG = cat(2,stackedEMG{:});
    
    

    disp('concatenating experiments...');
    CFA = [CFA; stackedCFA];
    RFA = [RFA; stackedRFA];
    full_EMG = cat(3,full_EMG,stackedEMG);
end
EMGtrials = mean(full_EMG,3);
disp('done!');

%% plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coef,score,latent] = pca(EMGtrials');
muscles_pc1 = score(:,1);

% remove zero rows
CFA( ~any(CFA,2), : ) = [];
RFA( ~any(RFA,2), : ) = [];

% remove low firing neurons
fr_cfa = sum(CFA,2)/(size(CFA,2)/1000);
fr_rfa = sum(RFA,2)/(size(RFA,2)/1000);
ex_cfa = find(fr_cfa<1);
ex_rfa = find(fr_rfa<1);

CFA(ex_cfa,:) = [];
RFA(ex_rfa,:) = [];


% normalize
normalize = @(x) (x-min(x))/max(x-min(x));

nCFA = CFA*0;
nRFA = RFA*0;
for i=1:size(nCFA,1)
    nCFA(i,:) = normalize(CFA(i,:));
end
for i=1:size(nRFA,1)
    nRFA(i,:) = normalize(RFA(i,:));
end

% get PCs of the neural activity
[~,score_CFA,latent_CFA] = pca(nCFA');
[~,score_RFA,latent_RFA] = pca(nRFA');

pc1_cfa = score_CFA(:,1);
pc2_cfa = score_CFA(:,2);
pc1_rfa = score_RFA(:,1);
pc2_rfa = score_RFA(:,2);


%
figure;
for i=1:3
    for j=1:8
        ind = 8*(i-1) + j;
        window = (200*(j-1))+1:j*200;
        subplot(3,8,ind)
%         disp([i j ind]);
        if(ind<=8)
            hold on
            plot(pc1_cfa(window));
            plot(pc2_cfa(window));
            hold off
            if(j==1)
                ylabel('first two PCs of CFA');
            end
            if(mod(j,2)==1)
                title(sprintf('reach spout %d', (j+1)/2));
            else
                title(sprintf('grasp spout %d', j/2));
            end
            xticks([1 100 200]);
            xticklabels([-100 0 100]);
        end
        if(ind>8 && ind<=16)
            hold on
            plot(pc1_rfa(window));
            plot(pc2_rfa(window));
            hold off
            if(j==1)
                ylabel('first two PCs of RFA');
            end
            xticks([1 100 200]);
            xticklabels([-100 0 100]);
        end
        if(ind>16 && ind<=24)
            plot(-100:99,muscles_pc1(window));    
            xline(0,':');
            if(j==1)
                ylabel('first principal component of EMG');
            end
        end
        
    end

end


disp(animal);















%%