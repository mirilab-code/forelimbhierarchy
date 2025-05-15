%% make spout trial plots



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
this_exp = exps{1};

totalCells1 = 0;
totalCells2 = 0;

disp(this_exp);
disp('loading things...');
datapath = fullfile('Z:\akiko\',this_exp,'\preprocess_with_acg');
files = fullfile(datapath,'*.mat');
matFiles = dir(files);
for i = 1:length(matFiles)
    baseFileName = fullfile(datapath, matFiles(i).name);
    load(baseFileName);
end

disp('getting trials bounds...');
statsIdx = aidx(1);
spoutNs = [stats.Spout1Used(statsIdx) stats.Spout2Used(statsIdx) stats.Spout3Used(statsIdx) stats.Spout4Used(statsIdx)];
disp(spoutNs);
spoutIdxs = cumsum(spoutNs);

spout_bounds = {};
spout_bounds{1} = reach_bounds_edit(1:spoutIdxs(1),:);
spout_bounds{2} = reach_bounds_edit(spoutIdxs(1)+1:spoutIdxs(2),:);
spout_bounds{3} = reach_bounds_edit(spoutIdxs(2)+1:spoutIdxs(3),:);
spout_bounds{4} = reach_bounds_edit(spoutIdxs(3)+1:spoutIdxs(4),:);


disp('done!');

%% rasters and behavior

% load cortical units from Adam
load('cortical_units/03272020_CFA_units.mat');
load('cortical_units/03272020_RFA_units.mat');
excludeCFA = ~ismember(events{1}(:,1),CFA_units);
events{1}(excludeCFA,:) = [];
excludeRFA = ~ismember(events{2}(:,1),RFA_units);
events{2}(excludeRFA,:) = [];

% change the indices so they don't skip a lot of numbers in the raster
units_cfa = unique(events{1}(:,1));
nUnits_cfa = length(units_cfa);
cfa_map = containers.Map(units_cfa,1:nUnits_cfa);
for i=1:length(events{1}(:,1))
    u = events{1}(i,1);
    events{1}(i,1) = cfa_map(u);
end
units_rfa = unique(events{2}(:,1));
nUnits_rfa = length(units_rfa);
rfa_map = containers.Map(units_rfa,1:nUnits_rfa);
for i=1:length(events{2}(:,1))
    u = events{2}(i,1);
    events{2}(i,1) = rfa_map(u);
end


muscles = {'Bicep','Tricep','ECR','PL'};
nMuscles = length(muscles);
all_reaches = sortrows(cat(1,spout_bounds{:}),1);

d = diff(all_reaches(:,1));

% look for 2 reaches that are close ish together and also are fairly
% fast/representative of a good reach

% reach1_ind = find(d==min(d));

%%
reach1_ind = 139    % a048 -> 24, 78 is good, 112 is good, 118 is good, 139 is MONEY 
reach2_ind = reach1_ind+1;

maxN = max(nUnits_rfa,nUnits_cfa);

% 

% reach1_ind = reach1_ind +1
% reach2_ind = reach2_ind +1

r1_start = all_reaches(reach1_ind,1);
r2_start = all_reaches(reach2_ind,1);
r1_end = all_reaches(reach1_ind,2);
r2_end = all_reaches(reach2_ind,2);

pad = 1000;

window = r1_start-pad:r2_end+(5*pad);
reach1_x = window(pad);
grasp1_x = window(pad+(r1_end-r1_start));
reach2_x = window(pad+r2_start-r1_start);
grasp2_x = window(pad+(r2_end-r1_start));


M = EMG(:,window);
for i=1:nMuscles
    z = zscore(M(i,:));
    M(i,:) = z;
end

%
clf('reset');
figure(1);
subplot(6,1,1);
hold on
plot_raster_from_events(events{1},2,window);
r = xline(reach1_x,'-r');
g = xline(grasp1_x,'-b');
xline(reach2_x,'-r');
xline(grasp2_x,'-b');
hold off
title('CFA');
axis tight
legend([r g],{'reach','grasp'});
ylim([0 maxN]);


subplot(6,1,2);
hold on
plot_raster_from_events(events{2},2,window);
xline(reach1_x,'-r');
xline(grasp1_x,'-b');
xline(reach2_x,'-r');
xline(grasp2_x,'-b');
ax = gca;
hold off
axis tight
title('RFA');
ylim([0 maxN]);


for i=1:nMuscles
    subplot(6,1,i+2);
    plot(M(i,:));
    ylabel(muscles{i});
    %     xline(reach1_x,'-r');
    %     xline(grasp1_x,'-b');
    %     xline(reach2_x,'-r');
    %     xline(grasp2_x,'-b');
    axis tight
end


%%