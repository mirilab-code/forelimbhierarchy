%% make the data necessary for CCA and PLS.

addpath(genpath('Z:/Scripts'));
base_folder = 'Z:\akiko\';

load('all_CFA_RFA_units.mat');

%%

folders = dir('Z:\akiko');
datafolders = struct2cell(folders);
foldernames = datafolders(1,3:end-10)';

opts = detectImportOptions('trialstats.csv');
opts = setvartype(opts, 'Date', 'char');
stats = readtable('trialstats.csv', opts);
dates = stats.Date(1:end-2);
animals = stats.Animal(1:end-2);


%% get just the cortical units for each session
animal_names = unique(animals,'stable');
all_exps = {};
for i=1:length(animal_names)
    animal = animal_names{i};
    aidx = find(contains(animals,animal));

    fidx = find(contains(foldernames,animal));
    all_exps = cat(1,all_exps,foldernames(fidx));
end

nexps = length(all_exps);

CFAunits = {};
RFAunits = {};
for i=1:nexps
    fname = [base_folder all_exps{i} '\preprocess_with_acg'];
    units_cfa = all_CFA_RFA_units{i}{1};
    units_rfa = all_CFA_RFA_units{i}{2};

    CFAunits{i} = units_cfa;
    RFAunits{i} = units_rfa;

end

%%

window = -99:100;
% lags = -60:60;
% lags = -100:100;
lags = [-500:10:-60 -50:50 60:10:500];

counter = 1;
for a=1:length(animal_names)
%%
%     if(a==2)    % note: you can't start from a certain number unless you mess with the counter. so just do all starting from 1.
%         break;
%     end

    animal = animal_names{a};
    aidx = find(contains(animals,animal));
    fidx = find(contains(foldernames,animal));
    exps = foldernames(fidx);
    for d=1:length(exps)

%         if(d==2)
%             break;
%         end

        disp(exps{d});
        disp('loading things...');
        datapath = fullfile('Z:\akiko\',exps{d},'\preprocess_with_acg');
        files = fullfile(datapath,'*.mat');
        matFiles = dir(files);
        for i = 1:length(matFiles)
            baseFileName = fullfile(datapath, matFiles(i).name);
            load(baseFileName);
        end

        % prune events to just have the cortical units
        cort_cfa = ismember(events{1}(:,1),CFAunits{counter});
        cort_rfa = ismember(events{2}(:,1),RFAunits{counter});
        events{1}(~cort_cfa,:) = [];
        events{2}(~cort_rfa,:) = [];

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
            trains{i} = trains{i}(~cellfun('isempty',trains{i}));
            FR{i} = train_to_firingrate(trains{i});
            FR{i} = FR{i}(:,1:duration);
        end

        % shift the spout_bounds by a random large number for the fake data
        fake_spout_bounds = spout_bounds;
        % not sure whether it's better to shift all the spout bounds by the
        % same number, or each a different number. Probably doesn't matter.
        % But technically it's possible to have the new bounds overlap if
        % they are each shifted by a different amount.
        % fake_spout_bounds = cellfun(@(x) x+randi([5*1000 10*1000],size(x,1),1), fake_spout_bounds, 'UniformOutput', false);
        fake_spout_bounds = cellfun(@(x) x+randi([5*1000 10*1000],1,1), fake_spout_bounds, 'UniformOutput', false);
        RFA_is_fake = true;

        CFA_RG = all_spouts_RG_lags(FR{1},spout_bounds,window,lags);
        RFA_RG = all_spouts_RG(FR{2},fake_spout_bounds,window);

        notnancolumns = logical(mean(~isnan(RFA_RG)));
        RFA_RG = RFA_RG(:,all(~isnan(RFA_RG)));   % for nan - columns
        
        CFA_RG = CFA_RG(:,notnancolumns,:);



        x = split(exps{d},'_');

        if(strcmp(x{1},'09232020'))
            disp('wait!')
        end

        savename = sprintf('output/%s_%s_CFAlagsandRFA_CONTROL', x{1},x{2});
        save(savename,'CFA_RG', 'RFA_RG','lags','RFA_is_fake');
        fprintf('%s: saved! \n', savename)

        counter = counter+1;
    end
end

% notes:
% the CFA_RG is a 3D matrix with depth being the lags. so for lags=-60:60,
% the zero lag index is 61. or just do lag0 = find(lags==0);

























%%