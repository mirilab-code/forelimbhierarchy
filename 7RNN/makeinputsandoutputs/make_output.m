%% 
addpath(genpath('Z:\code'));

%%
reaching_folder = 'Z:\reaching';
miceRFA = {'re7*','re8*','re9*'};
% miceRFA = {'re7*'};
miceCFA = {'re10*','re12*','re14*'};
% miceCFA = {'re10*'};

allmiceCFA = {};
for i=1:length(miceCFA)
    grp = get_all_subfolders(reaching_folder,miceCFA{i});
    allmiceCFA = cat(1,allmiceCFA,grp);
end
allmiceRFA = {};
for i=1:length(miceRFA)
    grp = get_all_subfolders(reaching_folder,miceRFA{i});
    allmiceRFA = cat(1,allmiceRFA,grp);
end

%%

% CFA inact mice
CFAinact = get_all_CL_avgs(allmiceCFA);
% RFA inact mice
RFAinact = get_all_CL_avgs(allmiceRFA);


disp('done!')
%%
disp('plotting');

window = RFAinact.window;

figure;
subplot(2,1,1);
hold on;
boundedline(window,CFAinact.MEAN_SUM_LASER,CFAinact.SEM_SUM_LASER,'c')
boundedline(window,CFAinact.MEAN_SUM_CONTROL,CFAinact.SEM_SUM_CONTROL,'k')
hold off
title('inact CFA')

subplot(2,1,2);
hold on;
boundedline(window,RFAinact.MEAN_SUM_LASER,RFAinact.SEM_SUM_LASER,'c')
boundedline(window,RFAinact.MEAN_SUM_CONTROL,RFAinact.SEM_SUM_CONTROL,'k')
hold off
title('inact RFA')
sgtitle('sum of muscle activity')




%%
muscle_names = {'bicep','tricep','ECR','PL'};
figure;
for i=1:4
    subplot(2,4,i);
    hold on;
    boundedline(window,CFAinact.MEAN_ALL_LASER(i,:),CFAinact.SEM_ALL_LASER(i,:),'c')
    boundedline(window,CFAinact.MEAN_ALL_CONTROL(i,:),CFAinact.SEM_ALL_CONTROL(i,:),'k')
    hold off
    title(muscle_names{i})
    if(i==1)
        ylabel('inactivate CFA')
    end
end


for i=1:4
    subplot(2,4,i+4);
    hold on;
    boundedline(window,RFAinact.MEAN_ALL_LASER(i,:),RFAinact.SEM_ALL_LASER(i,:),'c')
    boundedline(window,RFAinact.MEAN_ALL_CONTROL(i,:),RFAinact.SEM_ALL_CONTROL(i,:),'k')
    hold off
    title(muscle_names{i})
    if(i==1)
        ylabel('inactivate RFA')
    end
end



%% normalize or baseline subtract or whatever
bsln = 1:50;
rfa_control = RFAinact.MEAN_SUM_CONTROL - mean(RFAinact.MEAN_SUM_CONTROL(bsln));
cfa_control = CFAinact.MEAN_SUM_CONTROL - mean(CFAinact.MEAN_SUM_CONTROL(bsln));
rfa_laser = RFAinact.MEAN_SUM_LASER - mean(RFAinact.MEAN_SUM_LASER(bsln));
cfa_laser = CFAinact.MEAN_SUM_LASER - mean(CFAinact.MEAN_SUM_LASER(bsln));

control = mean([rfa_control; cfa_control]);

muscles.control = control;
muscles.inact_rfa = rfa_laser;
muscles.inact_cfa = cfa_laser;

figure;
hold on
plot(muscles.control)
plot(muscles.inact_rfa)
plot(muscles.inact_cfa)
hold off

%%


%% now do the neurons
load('model_setup_data\CFAinact.mat');
load('model_setup_data\RFAinact.mat');
neurons.no_stim_rfa = RFAinact.control_avg;
neurons.no_stim_cfa = CFAinact.control_avg;
neurons.inactRFA = RFAinact.laser_avg;
neurons.inactCFA = CFAinact.laser_avg;

neurons.no_stim_rfa = neurons.no_stim_rfa - mean(neurons.no_stim_rfa(bsln));
neurons.no_stim_cfa = neurons.no_stim_cfa - mean(neurons.no_stim_cfa(bsln));
neurons.inactRFA = neurons.inactRFA - mean(neurons.inactRFA(bsln));
neurons.inactCFA = neurons.inactCFA - mean(neurons.inactCFA(bsln));


hold on
plot(neurons.no_stim_cfa,'-r');
plot(neurons.no_stim_rfa,'-b');
plot(neurons.inactRFA,'-r');
plot(neurons.inactCFA,'-b');

%% 
output.muscles = muscles;
output.neurons = neurons;

save('model_setup_data\output.mat','output');


disp('done saving!');






%%
function X = get_all_CL_avgs(folder_list)


SUM_LASER = [];
ALL_LASER = [];
X.SEM_SUM_LASER  = [];
SUM_CONTROL = [];
ALL_CONTROL = [];
X.SEM_SUM_CONTROL = [];

for m=1:length(folder_list)
    data_folder = [folder_list{m} '\preprocess\'];
    files = fullfile(data_folder,'*.mat');
    matFiles = dir(files);
    for i = 1:length(matFiles)
        baseFileName = fullfile(data_folder, matFiles(i).name);
        load(baseFileName);
    end
    disp(data_folder);
    
    % z score the EMG
    EMG = zscore(EMG,0,2);

    laser = analogin(2,:);
    control = analogin(7,:);
    laser(laser<2) = 0;
    laser(laser>=2) = 1;
    control(control<2) = 0;
    control(control>=2) = 1;
    laser = logical(laser);
    control = logical(control);

    laser_on = find(diff(laser)>0);
    control_on = find(diff(control)>0);

    %%
    window = -50:100;

    L = get_trials(sum(EMG,1),laser_on,window);
    C = get_trials(sum(EMG,1),control_on,window);
    L_ALL = get_trials(EMG,laser_on,window);
    C_ALL = get_trials(EMG,control_on,window);
    

    SUM_LASER = cat(3,SUM_LASER,L);
    SUM_CONTROL = cat(3,SUM_CONTROL,C);
    ALL_LASER = cat(3,ALL_LASER,L_ALL);
    ALL_CONTROL = cat(3,ALL_CONTROL,C_ALL);


end

X.MEAN_SUM_LASER = mean(SUM_LASER,3);
X.SEM_SUM_LASER = std(SUM_LASER,0,3)/sqrt(size(SUM_LASER,3));
X.MEAN_SUM_CONTROL = mean(SUM_CONTROL,3);
X.SEM_SUM_CONTROL= std(SUM_CONTROL,0,3)/sqrt(size(SUM_CONTROL,3));

X.MEAN_ALL_LASER = mean(ALL_LASER,3);
X.SEM_ALL_LASER = std(ALL_LASER,0,3)/sqrt(size(ALL_LASER,3));
X.MEAN_ALL_CONTROL = mean(ALL_CONTROL,3);
X.SEM_ALL_CONTROL= std(ALL_CONTROL,0,3)/sqrt(size(ALL_CONTROL,3));




X.window = window;

end