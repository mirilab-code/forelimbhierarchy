%% get depth od probe histograms
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

load('cortical_units/03272020_CFA_units.mat');
load('cortical_units/03272020_RFA_units.mat');
excludeCFA = ~ismember(events{1}(:,1),CFA_units);
events{1}(excludeCFA,:) = [];
excludeRFA = ~ismember(events{2}(:,1),RFA_units);
events{2}(excludeRFA,:) = [];

disp('done!');

%%
depths_CFA = events{1}(:,3);
depths_CFA = sort(unique(depths_CFA));
depths_RFA = events{2}(:,3);
depths_RFA = sort(unique(depths_RFA));

depths_CFA = max(depths_CFA)-depths_CFA;
depths_RFA = max(depths_RFA)-depths_RFA;
% adjust for the angle 
depths_CFA = cosd(30)*depths_CFA;



mind = min(min(depths_CFA),min(depths_RFA));
maxd = max(max(depths_CFA),max(depths_RFA));
%%
bw = 100;

figure;
subplot(1,2,1)
histogram(depths_CFA,'BinWidth',bw);
xlabel('depth from pia (um)')
set(gca,'XDir','reverse');
camroll(90)
xlim([mind maxd*1.1]);
title('CFA')

subplot(1,2,2)
histogram(depths_RFA,'BinWidth',bw);
% set(gca,'XDir','reverse');
camroll(-90)
xlim([mind maxd*1.1]);
title('RFA')















%%