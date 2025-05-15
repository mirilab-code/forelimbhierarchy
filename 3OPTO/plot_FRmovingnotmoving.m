%% plot firing rates of neurons during movement and non movemet

addpath(genpath('Z:/code'))


base_folder = 'Z:\akiko\';


sessions = {...
    '03272020_a048_g0',...
    '03302020_a048_g0',...
    '03312020_a048_g0',...
    '04012020_a048_g0',...
    '04022020_a048_g0',...
    '04032020_a048_g0',...
    '09222020_a050_g0',...
    '09232020_a050_g0',...
    '09302020_ss2_g0',...
    '10012020_ss2_g0',...
    '10022020_ss2_g0',...
    '10052020_ss2_g0',...
    '10062020_ss2_g0',...
    '10242020_MA1_g0',...
    '10262020_MA1_g0',...
    '10272020_MA1_g0',...
    '10282020_MA1_g0',...
    '11102020_a051_g0',...
    '11122020_a051_g0',...
    '11152020_MA2_g0',...
    '11172020_MA2_g0',...
    };


nsessions = length(sessions);
%%

CELLS = struct;

for i=1:nsessions
    fname = ['muscleneurocorrelations\mayberight\' sessions{i} '.mat']
    load(fname);
    fname = ['muscleneurocorrelations\mayberight\' sessions{i} '_widthsanddepths.mat'];
    load(fname);
    events = load([base_folder sessions{i} '\preprocess_with_acg\events.mat']);
    events = events.events;
    duration = size(EMG,2);

    train_CFA = events_to_train(events{1});
    train_CFA = train_CFA(~cellfun('isempty',train_CFA));
    train_RFA = events_to_train(events{2});
    train_RFA = train_RFA(~cellfun('isempty',train_RFA));
    FR_cfa = trains_to_firingrate(train_CFA,duration);
    FR_rfa = trains_to_firingrate(train_RFA,duration);

    semg = sum(EMG,1);

    move = get_movement_from_RGbounds([base_folder sessions{i} '\preprocess_with_acg\'],semg);
    nomove = ~move;

    if(i==21)
        move = move(1:size(FR_cfa,2));
        nomove = ~move;
    end

    fr_move_CFA = FR_cfa(:,move);
    fr_nomove_CFA = FR_cfa(:,nomove);
    fr_move_RFA = FR_rfa(:,move);
    fr_nomove_RFA = FR_rfa(:,nomove);
   
    meanfr_move_CFA = mean(fr_move_CFA,2) * 1000;       % *1000 to make it in Hz
    meanfr_nomove_CFA = mean(fr_nomove_CFA,2) * 1000;
    meanfr_move_RFA = mean(fr_move_RFA,2) * 1000;
    meanfr_nomove_RFA = mean(fr_nomove_RFA,2) * 1000;

    CELLS(i).df = sessions{i};
    CELLS(i).meanfr_move_CFA = meanfr_move_CFA';
    CELLS(i).meanfr_nomove_CFA = meanfr_nomove_CFA';
    CELLS(i).wide_CFA = wide_CFA';
    CELLS(i).narrow_CFA = narrow_CFA';
    CELLS(i).meanfr_move_RFA = meanfr_move_RFA';
    CELLS(i).meanfr_nomove_RFA = meanfr_nomove_RFA';
    CELLS(i).wide_RFA = wide_RFA';
    CELLS(i).narrow_RFA = narrow_RFA';
end
save('movementnonmovement.mat',"CELLS");
%%
all_move_CFA = [CELLS.meanfr_move_CFA];
all_nomove_CFA = [CELLS.meanfr_nomove_CFA];
all_move_RFA = [CELLS.meanfr_move_RFA];
all_nomove_RFA = [CELLS.meanfr_nomove_RFA];
all_wide_CFA = [CELLS.wide_CFA];
all_narrow_CFA = [CELLS.narrow_CFA];
all_wide_RFA = [CELLS.wide_RFA];
all_narrow_RFA = [CELLS.narrow_RFA];

all_over1Hz_CFA = all_move_CFA>=1 | all_nomove_CFA>=1;
all_over1Hz_RFA = all_move_RFA>=1 | all_nomove_RFA>=1;

all_move_CFA_wide = all_move_CFA(all_wide_CFA);
all_nomove_CFA_wide = all_nomove_CFA(all_wide_CFA);
all_move_CFA_narrow = all_move_CFA(all_narrow_CFA);
all_nomove_CFA_narrow = all_nomove_CFA(all_narrow_CFA);

all_move_RFA_wide = all_move_RFA(all_wide_RFA);
all_nomove_RFA_wide = all_nomove_RFA(all_wide_RFA);
all_move_RFA_narrow = all_move_RFA(all_narrow_RFA);
all_nomove_RFA_narrow = all_nomove_RFA(all_narrow_RFA);

%%
lb = 0.1;
% total_max = max([all_move_CFA_wide all_move_CFA_narrow all_move_RFA_wide all_move_RFA_narrow])*1.5;
total_max = 100;

figure;
subplot(2,1,1);
hold on
scatter(all_move_CFA_wide,all_nomove_CFA_wide,'.r')
scatter(all_move_CFA_narrow,all_nomove_CFA_narrow,'.b')
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log')
refline(1,0)
xlim([lb total_max])
ylim([lb total_max])
title('CFA')
legend({'wide','narrow'},'Location','northwest')
set(gca, TickDir='out');

subplot(2,1,2);
hold on
scatter(all_move_RFA_wide,all_nomove_RFA_wide,'.r')
scatter(all_move_RFA_narrow,all_nomove_RFA_narrow,'.b')
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log')
refline(1,0)
xlim([lb total_max])
ylim([lb total_max])
title('RFA')

set(gca, TickDir='out');

xlabel('firing rate during movement')
ylabel('firing rate at rest')


%%




%%
function movement = get_movement_from_RGbounds(df,vec)
    files = dir(df);
    names = {files.name};
    ind = cellfun(@(x) contains(x,'reach_bounds_edit'),names);
%     disp(names{ind});
    
    window = -99:100;
    RGbounds = load([df names{ind}]);
    RGbounds = RGbounds.reach_bounds_edit;

    movement = logical(vec*0);
    reach = window + RGbounds(:,1);
    reach = reach(:);
    grasp = window + RGbounds(:,2);
    grasp = grasp(:);

    movement(reach) = 1;
    movement(grasp) = 1;




end





















%%
