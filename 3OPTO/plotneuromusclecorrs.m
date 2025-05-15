%% find and plot neurons that correlate with muscles


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


%% maybe trash
CELLS = struct;

for i=1:nsessions
    fname = ['muscleneurocorrelations\mayberight\' sessions{i} '.mat']
    load(fname);
    fname = ['muscleneurocorrelations\mayberight\' sessions{i} '_widthsanddepths.mat'];
    load(fname);
    CELLS(i).datafolder = sessions{i};

    semg = sum(EMG,1);
    duration = size(EMG,2);

    % firing rate for each region
    fr_CFA = cellfun(@(x) length(x)/(duration/1000),train_CFA);
    fr_RFA = cellfun(@(x) length(x)/(duration/1000),train_RFA);
    CELLS(i).firingrateCFA = fr_CFA';
    CELLS(i).firingrateRFA = fr_RFA';
    over1Hz_CFA = [fr_CFA>=1]';
    over1Hz_RFA = [fr_RFA>=1]';
    CELLS(i).over1Hz_CFA = over1Hz_CFA;
    CELLS(i).over1Hz_RFA = over1Hz_RFA;


    % depths
    CELLS(i).depthsCFA = depths_CFA;
    CELLS(i).depthsRFA = depths_RFA;

    % widths
    CELLS(i).narrowCFA = narrow_CFA;
    CELLS(i).narrowRFA = narrow_RFA;
    CELLS(i).wideCFA = wide_CFA;
    CELLS(i).wideRFA = wide_RFA;

    % significantly correlated neurons
    frac_sig_CFA = mean(sigCFA);
    frac_sig_RFA = mean(sigRFA);
    frac_sig_CFA_over1Hz = mean(sigCFA(over1Hz_CFA));
    frac_sig_RFA_over1Hz = mean(sigRFA(over1Hz_RFA));
    CELLS(i).fracsigCFA = frac_sig_CFA;
    CELLS(i).fracsigRFA = frac_sig_RFA;
    CELLS(i).fracsigCFA_over1Hz = frac_sig_CFA_over1Hz;
    CELLS(i).fracsigRFA_over1Hz = frac_sig_RFA_over1Hz;



    clearvars -except sessions nsessions CELLS i
end



%% plotting


% frac_sig = [[CELLS.fracsigCFA]' [CELLS.fracsigRFA]'];
frac_sig = [[CELLS.fracsigCFA_over1Hz]' [CELLS.fracsigRFA_over1Hz]'];
figure;
title('fraction of neurons that correlate with movement')
hold on;
scatter(frac_sig(:,1),frac_sig(:,2),'ok','filled');
ylim([0 1]);
xlim([0 1]);
rl = refline(1,0);
rl.Color = 'black';
rl.LineStyle = '--';
xlabel('CFA')
ylabel('RFA')

%%


figure;
hold on
plot([0 1],frac_sig','-k')
scatter(frac_sig(:,1)*0,frac_sig(:,1),'or','filled')

scatter(0,mean(frac_sig(:,1)),'dk','filled')

scatter(frac_sig(:,2)*0+1,frac_sig(:,2),'ob','filled')
scatter(1,mean(frac_sig(:,2)),'dk','filled')
xlim([-1 2])

xticklabels({'','','CFA','','RFA',''})
ylabel('fraction of neurons')
title('significcantly correlated with at least one muscle and fire over 1Hz')
set(gca, "TickDir",'out');

%% now do it but aggregate across animals weighting by the number of recorded neurons
nCFA = {CELLS.depthsCFA};
nCFA = cellfun(@length,nCFA);
nRFA = {CELLS.depthsRFA};
nRFA = cellfun(@length,nRFA);

a048 = [1 2 3 4 5 6];
a050 = [7 8];
ss2 = [10 11 12 13];
MA1 = [14 15 16 17];
a051 = [18 19];
MA2 = [20 21];

frac_sig_over1Hz_CFA = [CELLS.fracsigCFA_over1Hz];
frac_sig_over1Hz_RFA = [CELLS.fracsigRFA_over1Hz];

mice = {a048,a050,ss2,MA1,a051,MA2};
frac_sig_weighted = [];
for i=1:length(mice)
    this_mouse = mice{i};
    this_CFA = nCFA(this_mouse);
    this_RFA = nRFA(this_mouse);
    this_FS_CFA = frac_sig_over1Hz_CFA(this_mouse);
    this_FS_RFA = frac_sig_over1Hz_RFA(this_mouse);
    wCFA = this_CFA / sum(this_CFA);
    wRFA = this_RFA / sum(this_RFA);


    FSW_CFA = dot(wCFA,this_FS_CFA);
    FSW_RFA = dot(wRFA,this_FS_RFA);
    frac_sig_weighted = [frac_sig_weighted; [FSW_CFA FSW_RFA]];

end
disp('done!')
%
figure;
hold on
plot([0 1],frac_sig_weighted','-k')
scatter(frac_sig_weighted(:,1)*0,frac_sig_weighted(:,1),'or','filled')

scatter(0,mean(frac_sig_weighted(:,1)),'dk','filled')

scatter(frac_sig_weighted(:,2)*0+1,frac_sig_weighted(:,2),'ob','filled')
scatter(1,mean(frac_sig_weighted(:,2)),'dk','filled')
xlim([-1 2])

xticklabels({'','','CFA','','RFA',''})
ylabel('fraction of neurons')
title('significcantly correlated with at least one muscle and fire over 1Hz')
set(gca, "TickDir",'out');

%% firing rates by animal and total
% all firing rates

fr_CFA = {CELLS.firingrateCFA};
fr_CFA_narrow = {CELLS.narrowCFA};
fr_CFA_wide = {CELLS.wideCFA};
fr_CFA_total = [CELLS.firingrateCFA];


fr_RFA = {CELLS.firingrateRFA};
fr_RFA_narrow = {CELLS.narrowRFA};
fr_RFA_wide = {CELLS.wideRFA};
fr_RFA_total = [CELLS.firingrateRFA];

figure;
ax1 = subplot(2,2,1);
hold on
plot_firingratehist(ax1,fr_CFA,fr_CFA_wide,fr_CFA_total,fr_RFA_total);
title('CFA wide')

ax2 = subplot(2,2,2);
title('CFA narrow')
hold on
plot_firingratehist(ax2,fr_CFA,fr_CFA_narrow,fr_CFA_total,fr_RFA_total);

ax3 = subplot(2,2,3);
title('RFA wide')
hold on
plot_firingratehist(ax3,fr_RFA,fr_RFA_wide,fr_CFA_total,fr_RFA_total);

ax4 = subplot(2,2,4);
title('RFA narrow')
hold on
plot_firingratehist(ax4,fr_RFA,fr_RFA_narrow,fr_CFA_total,fr_RFA_total);

linkaxes([ax1 ax2],'xy')

%%
function plot_firingratehist(ax,fr_cell,width_type_cell,all_fr_cfa,all_fr_rfa)

maxfr = max([all_fr_cfa(:); all_fr_rfa(:)]);
edges = 0:0.1:maxfr;
all_fr = [fr_cell{:}];
all_widthtype = [];

for i=1:length(fr_cell)
    this_fr = fr_cell{i};
    wt = width_type_cell{i};
    all_widthtype = [all_widthtype; wt];
    h = histcounts(this_fr(wt),BinEdges=edges, Normalization='cdf');
    plot(edges(1:end-1),h,Color=[0 0 0 0.5], Parent=ax);
end
all_widthtype = logical(all_widthtype);
[h] = histcounts(all_fr(all_widthtype),BinEdges=edges, Normalization='cdf');
plot(edges(1:end-1),h,Color=[1 0 0 1]);
set(gca,'Xscale','log');
xlabel('log firing rate');
ylabel('fraction of units')

end

%%


%% this is the exact same result as above
% CELLS = struct;
% 
% for i=1:nsessions
%     fname = ['muscleneurocorrelations\mayberight\' sessions{i} '.mat']
%     load(fname);
%     fname = ['muscleneurocorrelations\mayberight\' sessions{i} '_widthsanddepths.mat'];
%     load(fname);
%     events = load([base_folder sessions{i} '\preprocess_with_acg\events.mat']);
%     events = events.events;
%     duration = size(EMG,2);
% 
%     train_CFA = events_to_train(events{1});
%     train_CFA = train_CFA(~cellfun('isempty',train_CFA));
%     train_RFA = events_to_train(events{2});
%     train_RFA = train_RFA(~cellfun('isempty',train_RFA));
%     FR_cfa = trains_to_firingrate(train_CFA,duration);
%     FR_rfa = trains_to_firingrate(train_RFA,duration);
% 
%     obs_corr_CFA = corr(FR_cfa',EMG');  
%     obs_corr_RFA = corr(FR_rfa',EMG');  
% 
%     bcfa = permute(bootstrap_CFA,[2 1 3]);
%     brfa = permute(bootstrap_RFA,[2 1 3]);
% 
%     % find the fraction of observed correlations that are greater/less than the bootstrapped null distributions
%     frac_CFA_top = mean(obs_corr_CFA > bcfa,3);
%     frac_CFA_bottom = mean(obs_corr_CFA < bcfa,3);
%     frac_RFA_top = mean(obs_corr_RFA > brfa,3);
%     frac_RFA_bottom = mean(obs_corr_RFA < brfa,3);
% 
%     % two tailed test so test if things are above or below 2.5% of the population
%     sig_CFA_top = frac_CFA_top > (1-0.0025);
%     sig_CFA_bottom = frac_CFA_bottom < 0.0025;
%     sig_CFA = sig_CFA_bottom | sig_CFA_top;
%     sig_CFA = any(sig_CFA,2);
%     sig_RFA_top = frac_RFA_top > (1-0.0025);
%     sig_RFA_bottom = frac_RFA_bottom < 0.0025;
%     sig_RFA = sig_RFA_bottom | sig_RFA_top;
%     sig_RFA = any(sig_RFA,2);
% 
%     fracsigCFA = mean(sig_CFA);
%     fracsigRFA = mean(sig_RFA);
% 
%     CELLS(i).datafolder = sessions{i};
% 
%     % depths
%     CELLS(i).depthsCFA = depths_CFA;
%     CELLS(i).depthsRFA = depths_RFA;
% 
%     % widths
%     CELLS(i).narrowCFA = narrow_CFA;
%     CELLS(i).narrowRFA = narrow_RFA;
%     CELLS(i).wideCFA = wide_CFA;
%     CELLS(i).wideRFA = wide_RFA;
% 
%     % significantly correlated neurons
%     frac_sig_CFA = mean(sigCFA);
%     frac_sig_RFA = mean(sigRFA);
%     CELLS(i).fracsigCFA = frac_sig_CFA;
%     CELLS(i).fracsigRFA = frac_sig_RFA;
% end




%%