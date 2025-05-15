%% make figure plots


addpath(genpath('Z:/code'));
disp('unconstrained first');
data_folder_unconstrained = uigetdir('C:\Users\mirilab\Documents\Mark\pythonRNN\output_unconstrained');
act_files = dir([data_folder_unconstrained '\activity_*.npy']);
As_unconstrained = {act_files.name};
network_files = dir([data_folder_unconstrained '\newW_*.npy']);
Ws_unconstrained = {network_files.name};
unconstrained = get_model_data(Ws_unconstrained,As_unconstrained,data_folder_unconstrained);

%
disp('then symmetric');
data_folder_symmetric = uigetdir('C:\Users\mirilab\Documents\Mark\pythonRNN\output_symmetric');
act_files = dir([data_folder_symmetric '\activity_*.npy']);
As_symmetric = {act_files.name};
network_files = dir([data_folder_symmetric '\newW_*.npy']);
Ws_symmetric = {network_files.name};
symmetric = get_model_data(Ws_symmetric,As_symmetric,data_folder_symmetric);



%% unconstrained interregional weight histograms
figure(1);
hold on; 
histogram(unconstrained.total_CtoR,'BinWidth',0.025,'Normalization','probability');
histogram(unconstrained.total_RtoC,'BinWidth',0.025,'Normalization','probability');

%% do it but log scale and also box plots
log_CtoR = log(unconstrained.total_CtoR);
log_CtoR(log_CtoR<-10) = [];
log_RtoC = log(unconstrained.total_RtoC);
log_RtoC(log_RtoC<-10) = [];

bw = 0.5;
figure;
hold on; 
histogram(log_CtoR,'BinWidth',bw,'Normalization','probability');
histogram(log_RtoC,'BinWidth',bw,'Normalization','probability');
legend('CtoR','RtoC')
title(bw)

bw = 0.25;
figure;
hold on; 
histogram(log_CtoR,'BinWidth',bw,'Normalization','probability');
histogram(log_RtoC,'BinWidth',bw,'Normalization','probability');
legend('CtoR','RtoC')
title(bw)

bw = 1;
figure;
hold on; 
histogram(log_CtoR,'BinWidth',bw,'Normalization','probability');
histogram(log_RtoC,'BinWidth',bw,'Normalization','probability');
legend('CtoR','RtoC')
title(bw)

%% then also do box plots
figure;
boxplot(log_CtoR,'Whisker',Inf)
title('CFA to RFA')
ylim([-11 3])
ylabel('log weights')

figure;
boxplot(log_RtoC,'Whisker',Inf)
title('RFA to CFA')
ylim([-11 3])
ylabel('log weights')
%% unconstrained weighted activity of interregional units
mean_unc_act_from_c = squeeze(mean(unconstrained.no_stim_from_c,2));
mean_unc_act_from_r = squeeze(mean(unconstrained.no_stim_from_r,2));

% -25 to 0, := onset of go cue (go cue = halfmax at t=-25 from reach onset) to reach onset
act_neg25to0_from_c = mean(mean_unc_act_from_c(25:50,:),1);
act_neg25to0_from_r = mean(mean_unc_act_from_r(25:50,:),1);

% 0 to 50, from reach onset
act_0to50_from_c = mean(mean_unc_act_from_c(50:100,:),1);
act_0to50_from_r = mean(mean_unc_act_from_r(50:100,:),1);

% 50 to 100, from reach onset
act_50to100_from_c = mean(mean_unc_act_from_c(100:150,:),1);
act_50to100_from_r = mean(mean_unc_act_from_r(100:150,:),1);

figure(2);
subplot(1,2,1);
hold on
boxplot([act_neg25to0_from_c; act_neg25to0_from_r; act_0to50_from_c; act_0to50_from_r; act_50to100_from_c; act_50to100_from_r]', ...
    'whisker',Inf,'Labels',{'-25 to 0 from CFA','-25 to 0 from RFA','0 to 50 from CFA','0 to 50 from RFA','50 to 100 from CFA','50 to 100 from RFA'})
ylabel('weighted activity')

subplot(1,2,2);
hold on
boxplot([act_neg25to0_from_c; act_0to50_from_c; act_50to100_from_c; act_neg25to0_from_r; act_0to50_from_r; act_50to100_from_r]', ...
    'whisker',Inf,'Labels', ...
    {'-25 to 0 from CFA','0 to 50 from CFA','50 to 100 from CFA','-25 to 0 from RFA','0 to 50 from RFA','50 to 100 from RFA'})
sgtitle('unconstrained weighted activty. these are the same but the order is different, idk which is better')
%% performance of both regimes
figure(3)
hold on
boxplot([unconstrained.performance symmetric.performance],'Whisker',Inf,'Labels',{'Unconstrained','Symmetric'})
hold off
ylabel('performance')
ylim([0 250])
title('performance of the two types of models')

%% symmetric: weighted activity of interregional units
mean_sym_act_from_c = squeeze(mean(symmetric.no_stim_from_c,2));
mean_sym_act_from_r = squeeze(mean(symmetric.no_stim_from_r,2));

% -25 to 0, := onset of go cue (go cue = halfmax at t=-25 from reach onset) to reach onset
act_neg25to0_from_c = mean(mean_sym_act_from_c(25:50,:),1);
act_neg25to0_from_r = mean(mean_sym_act_from_r(25:50,:),1);

% 0 to 50, from reach onset
act_0to50_from_c = mean(mean_sym_act_from_c(50:100,:),1);
act_0to50_from_r = mean(mean_sym_act_from_r(50:100,:),1);

% 50 to 100, from reach onset
act_50to100_from_c = mean(mean_sym_act_from_c(100:150,:),1);
act_50to100_from_r = mean(mean_sym_act_from_r(100:150,:),1);



figure(4);
subplot(1,2,1);
hold on
boxplot([act_neg25to0_from_c; act_neg25to0_from_r; act_0to50_from_c; act_0to50_from_r; act_50to100_from_c; act_50to100_from_r]', ...
    'whisker',Inf,'Labels',{'-25 to 0 from CFA','-25 to 0 from RFA','0 to 50 from CFA','0 to 50 from RFA','50 to 100 from CFA','50 to 100 from RFA'})
ylabel('weighted activity')

subplot(1,2,2);
hold on
boxplot([act_neg25to0_from_c; act_0to50_from_c; act_50to100_from_c; act_neg25to0_from_r; act_0to50_from_r; act_50to100_from_r]', ...
    'whisker',Inf,'Labels', ...
    {'-25 to 0 from CFA','0 to 50 from CFA','50 to 100 from CFA','-25 to 0 from RFA','0 to 50 from RFA','50 to 100 from RFA'})
sgtitle('symmetric weighted activty. these are the same but the order is different, idk which is better')

%% intre regional exc and inh weights

figure(5);
subplot(1,2,1);
hold on
boxplot([unconstrained.avg_CtoC_exc unconstrained.avg_RtoR_exc -unconstrained.avg_CtoC_inh -unconstrained.avg_RtoR_inh],'Whisker',Inf, ...
    'Labels',{'CFA exc','RFA exc','CFA inh','RFA inh'});
ylim([0 0.2])
title('unconstrained intraregional weights')

subplot(1,2,2);
hold on
boxplot([symmetric.avg_CtoC_exc symmetric.avg_RtoR_exc -symmetric.avg_CtoC_inh -symmetric.avg_RtoR_inh],'Whisker',Inf, ...
    'Labels',{'CFA exc','RFA exc','CFA inh','RFA inh'});
ylim([0 0.2])
title('symmetric intraregional weights')


%% do some t tests

% excitatory unconstrained CFA intra <-> RFA intra
[h_EU,p_EU] = ttest2(unconstrained.avg_CtoC_exc,unconstrained.avg_RtoR_exc);

% inhibitory unconstrained CFA intra <-> RFA intra
[h_IU,p_IU] = ttest2(unconstrained.avg_CtoC_inh,unconstrained.avg_RtoR_inh)

% excitatory constrained CFA intra <-> RFA intra
[h_EC,p_EC] = ttest2(symmetric.avg_CtoC_exc,symmetric.avg_RtoR_exc);

% inhibitory constrained CFA intra <-> RFA intra
[h_IC,p_IC] = ttest2(symmetric.avg_CtoC_inh,symmetric.avg_RtoR_inh)








%%
function X = get_model_data(Ws,As,data_folder)

W = zeros(1000,1000,length(Ws));

X.no_stim_from_c = [];
X.measurecfa_inactrfa_from_c = [];
X.measurerfa_inactcfa_from_c = [];
X.no_stim_from_r = [];
X.measurecfa_inactrfa_from_r = [];
X.measurerfa_inactcfa_from_r = [];
X.ctoc_sum = [];
X.rtor_sum = [];

X.total_CtoR = [];
X.total_RtoC = [];

X.sum_CtoR = [];
X.avg_CtoC_exc = [];
X.avg_CtoC_inh = [];
X.var_CtoC_exc = [];
X.var_CtoC_inh = [];
X.sum_RtoC = [];
X.avg_RtoR_exc = [];
X.avg_RtoR_inh = [];
X.var_RtoR_exc = [];
X.var_RtoR_inh = [];
X.Ws =[];

X.performance = [];
for i=1:length(As)
    act_file = As{i};
    w_file = Ws{i};
    perf = split(w_file,'_');
    perf = perf{end};
    perf = perf(1:end-4);   % get rid of the .npy
    perf = str2num(perf);
    X.performance = [X.performance; perf];
    act = readNPY([data_folder '\' act_file]);
    w = readNPY([data_folder '\' w_file]);
    X.Ws = cat(3,X.Ws,w);

    % intRAregional weights

    % first just sum together the intra weights
    X.ctoc_sum = [X.ctoc_sum; sum(sum(w(1:500,1:500)))];
    X.rtor_sum = [X.rtor_sum; sum(sum(w(501:end,501:end)))];

    ctoc_exc = w(1:500,1:400);
    ctoc_exc = ctoc_exc(:);
    ctoc_exc(ctoc_exc==0) = [];
    ctoc_inh = w(1:500,401:end);
    ctoc_inh = ctoc_inh(:);
    ctoc_inh(ctoc_inh==0) = [];
    
    rtor_exc = w(501:end,501:900);
    rtor_exc = rtor_exc(:);
    rtor_exc(rtor_exc==0) = [];
    rtor_inh = w(501:end,901:end);
    rtor_inh = rtor_inh(:);
    rtor_inh(rtor_inh==0) = [];

    X.avg_CtoC_exc = [X.avg_CtoC_exc; mean(ctoc_exc)];
    X.avg_RtoR_exc = [X.avg_RtoR_exc; mean(rtor_exc)];
    X.avg_CtoC_inh = [X.avg_CtoC_inh; mean(ctoc_inh)];
    X.avg_RtoR_inh = [X.avg_RtoR_inh; mean(rtor_inh)];

    X.var_CtoC_exc = [X.var_CtoC_exc; var(ctoc_exc)];
    X.var_RtoR_exc = [X.var_RtoR_exc; var(rtor_exc)];
    X.var_CtoC_inh = [X.var_CtoC_inh; var(ctoc_inh)];
    X.var_RtoR_inh = [X.var_RtoR_inh; var(rtor_inh)];

    CtoR_w = w(501:1000,1:500);
    RtoC_w = w(1:500,501:1000);
    X.sum_CtoR = [X.sum_CtoR; sum(CtoR_w(:))];
    X.sum_RtoC = [X.sum_RtoC; sum(RtoC_w(:))];
    X.total_CtoR = [X.total_CtoR; nonzeros(CtoR_w(:))];
    X.total_RtoC = [X.total_RtoC; nonzeros(RtoC_w(:))];

    inter_from_c = any(CtoR_w,1);
    inter_from_r = any(RtoC_w,1);

    projection_weights_from_c = sum(CtoR_w(:,inter_from_c));
    projection_weights_from_r = sum(RtoC_w(:,inter_from_r));

    no_stim_from_c = act(:,inter_from_c,1) .* projection_weights_from_c;
    measurecfa_inactrfa_from_c = act(:,inter_from_c,2) .* projection_weights_from_c;
    measurerfa_inactcfa_from_c = act(:,inter_from_c,3) .* projection_weights_from_c;
    no_stim_from_r = act(:,inter_from_r,1) .* projection_weights_from_r;
    measurecfa_inactrfa_from_r = act(:,inter_from_r,2) .* projection_weights_from_r;
    measurerfa_inactcfa_from_r = act(:,inter_from_r,3) .* projection_weights_from_r;

    X.no_stim_from_c = cat(3,X.no_stim_from_c,no_stim_from_c);
    X.measurecfa_inactrfa_from_c = cat(3,X.measurecfa_inactrfa_from_c,measurecfa_inactrfa_from_c);
    X.measurerfa_inactcfa_from_c = cat(3,X.measurerfa_inactcfa_from_c,measurerfa_inactcfa_from_c);

    X.no_stim_from_r = cat(3,X.no_stim_from_r,no_stim_from_r);
    X.measurecfa_inactrfa_from_r = cat(3,X.measurecfa_inactrfa_from_r,measurecfa_inactrfa_from_r);
    X.measurerfa_inactcfa_from_r = cat(3,X.measurerfa_inactcfa_from_r,measurerfa_inactcfa_from_r);

end
X.CminusR = X.sum_CtoR - X.sum_RtoC;  %positive values = more CtoR


end
