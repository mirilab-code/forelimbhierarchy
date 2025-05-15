function [sname] = getwidthsanddepths(data_folder)

% disp('loading things...');
disp(data_folder);
files = fullfile(data_folder,'*.mat');
matFiles = dir(files);


% do neurons
events_ind = find(strcmp({matFiles.name},'events.mat'));
events_fname = fullfile(data_folder, matFiles(events_ind).name);
data = load(events_fname);
events = data.events;
train_CFA = events_to_train(events{1});
train_RFA = events_to_train(events{2});
nonemptyCFA = ~cellfun('isempty',train_CFA);
nonemptyRFA = ~cellfun('isempty',train_RFA);

% widths
widths_ind = find(strcmp({matFiles.name},'widths.mat'));
widths_fname = fullfile(data_folder, matFiles(widths_ind).name);
widths = load(widths_fname);
widths_CFA  = widths.spike_widths{1};
widths_RFA  = widths.spike_widths{2};

widths_CFA = widths_CFA(1:length(nonemptyCFA)); % uhhh...
widths_RFA = widths_RFA(1:length(nonemptyRFA));

% then emg
emg_ind = find(strcmp({matFiles.name},'EMG.mat'));
emg_fname = fullfile(data_folder, matFiles(emg_ind).name);
data = load(emg_fname);

EMG = data.EMG;


train_CFA = train_CFA(nonemptyCFA);
train_RFA = train_RFA(nonemptyRFA);
widths_CFA = widths_CFA(nonemptyCFA);
widths_RFA = widths_RFA(nonemptyRFA);


disp('done loading')
%%
narrow_CFA = widths_CFA <= 13;
wide_CFA = ~narrow_CFA;
narrow_RFA = widths_RFA <= 13;
wide_RFA = ~narrow_RFA;

cfa_inds = events{1}(:,1);
cfa_depths = events{1}(:,3);
cfa_m = containers.Map(cfa_inds,cfa_depths);
rfa_inds = events{2}(:,1);
rfa_depths = events{2}(:,3);
rfa_m = containers.Map(rfa_inds,rfa_depths);

depths_CFA = values(cfa_m);
depths_RFA = values(rfa_m);









%%
sname = split(data_folder,'\');
sname = ['muscleneurocorrelations\' sname{end-1} '_widthsanddepths.mat'];

save(sname,'narrow_CFA','wide_CFA','narrow_RFA','wide_RFA','depths_CFA','depths_RFA','widths_CFA','widths_RFA');
fprintf('saved to %s \n',sname);




end