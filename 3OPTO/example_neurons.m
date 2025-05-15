%% get example neurons and plot pre spout
addpath(genpath('Z:\Scripts\'));

% load stuff
data_folder = 'Z:\akiko\10052020_ss2_g0\preprocess_with_acg';

idcs   = strfind(data_folder,'\');
files = fullfile(data_folder,'*.mat');
matFiles = dir(files);
for i = 1:length(matFiles)
    baseFileName = fullfile(data_folder, matFiles(i).name);
    load(baseFileName);
end

disp(data_folder);
CFA_events = events{1};
RFA_events = events{2};

duration = size(EMG,2);
CFA_train = events_to_train(CFA_events);
RFA_train = events_to_train(RFA_events);

fr_CFA = cellfun(@(x) length(x)/duration, CFA_train);
fr_RFA = cellfun(@(x) length(x)/duration, RFA_train);

%% make a janky neurons data structure with train|firing rate|width
CFA = struct;
ucfa = 1:length(CFA_train);
for i=1:length(ucfa)
    CFA(i).train = CFA_train{i};
    CFA(i).width = spike_widths{1}(i);
    CFA(i).firingrate = fr_CFA(i);
end
empty = cellfun(@isempty, CFA_train);
CFA(empty) = [];

RFA = struct;
urfa = 1:length(RFA_train);
for i=1:length(urfa)
    RFA(i).train = RFA_train{i};
    RFA(i).width = spike_widths{2}(i);
    RFA(i).firingrate = fr_RFA(i);
end
empty = cellfun(@isempty, RFA_train);
RFA(empty) = [];

%% find the highest firing rate neurons from each region. we want 2 pyramidal (wide) and 1 interneuron (narrow) from each.
cfa_fr = [CFA.firingrate];
cfa_width = [CFA.width];
rfa_fr = [RFA.firingrate];
rfa_width = [RFA.width];
%%
cfa_pyr1 = 171; % keep I guess
cfa_pyr2 = 212; % def keep
cfa_inter = 159;    %def keep
rfa_pyr1 = 33;      % for RFA: pyr could be 27,34,41,44      inter could be 27,32, 
rfa_pyr2 = 85;  % maybe keep 28         
rfa_inter = 42; %def keep 42


%% extract spike trains for reach and grasp at each spout
window = -200:200;

u = [cfa_pyr1; cfa_pyr2; cfa_inter; rfa_pyr1; rfa_pyr2; rfa_inter];
types = {'pyramidal CFA','pyramidal CFA','interneuron CFA','pyramidal RFA','pyramidal RFA','interneuron RFA'};

raster_trials = cell(6,1);
for n=1:3
    raster_trials{n} = cell(4,1);
    for i=1:4
        raster_trials{n}{i} = cell(2,1);
        raster_trials{n}{i}{1} = get_raster_trials(CFA(u(n)).train,reach_grasp{i}(:,1),window);
        raster_trials{n}{i}{2} = get_raster_trials(CFA(u(n)).train,reach_grasp{i}(:,2),window);
    end
end
for n=4:6
    raster_trials{n} = cell(4,1);
    for i=1:4
        raster_trials{n}{i} = cell(2,1);
        raster_trials{n}{i}{1} = get_raster_trials(RFA(u(n)).train,reach_grasp{i}(:,1),window);
        raster_trials{n}{i}{2} = get_raster_trials(RFA(u(n)).train,reach_grasp{i}(:,2),window);
    end
end
disp('done getting trials')

%%
x = {raster_trials{1}{:}};
x = vertcat(x{:});
n_trials_spouts = cellfun(@length,x);
mintrials = min(n_trials_spouts);

%%
mkrsz = 1;
for n=1:6
    f(n) = figure(n);
    tiledlayout(2,4)
    for i=1:2
        for j=1:4
            
            nexttile;

            T = raster_trials{n}{j}{i};
            T = T(randperm(numel(T)));
            T = T(1:mintrials);
            hold on
            for trn=1:length(T)
                x = ones(1,length(T{trn}))*trn;
                scatter(T{trn},x, mkrsz, 'k', '|');
            end
            hold off
            axis tight
            ylim([0 length(T)]);

            if(i==1 && j==1)
                ylabel('reach')
            elseif(i==2 && j==1)
                ylabel('grasp')
            end
            if(i==1)
                ss = sprintf('spout %d',j);
                title(ss)
            end
            xlim([window(1) window(end)])
        end
    end
    sgtitle(sprintf('unit %d %s',u(n),types{n}));
    fname = sprintf('unit %d %s',u(n),types{n});
    print('-painters','-dsvg',fname)
end

%%
ylims = [120 160 80 50 80 60];
bw = 20;
for n=1:6
    f(n) = figure(n);
    tiledlayout(2,4)
    for i=1:2
        for j=1:4
            
            nexttile;

            T = raster_trials{n}{j}{i};
            T = T(randperm(numel(T)));
            T = T(1:mintrials);
            hold on
            h = histogram(cell2mat(T),'BinWidth',bw);
            hold off
            axis tight

            ylim([0 ylims(n)]);


            if(i==1 && j==1)
                ylabel('reach')
            elseif(i==2 && j==1)
                ylabel('grasp')
            end
            if(i==1)
                ss = sprintf('spout %d',j);
                title(ss)
            end
            xlim([window(1) window(end)])
        end
    end
    sgtitle(sprintf('unit %d %s',u(n),types{n}));
    fname = sprintf('FR histogram unit %d %s',u(n),types{n});
    print('-painters','-dsvg',fname)
end

















%%