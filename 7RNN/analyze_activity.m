%% analyze activity


addpath(genpath('Z:/code'));
disp('unconstrained first');
data_folder_unconstrained = uigetdir('C:\Users\mirilab\Documents\Mark\pythonRNN\output_unconstrained');
act_files = dir([data_folder_unconstrained '\activity_*.npy']);
As_unconstrained = {act_files.name};
network_files = dir([data_folder_unconstrained '\newW_*.npy']);
Ws_unconstrained = {network_files.name};
unconstrained = get_model_data(Ws_unconstrained,As_unconstrained,data_folder_unconstrained);

%%
disp('then symmetric');
data_folder_symmetric = uigetdir('C:\Users\mirilab\Documents\Mark\pythonRNN\output_symmetric');
act_files = dir([data_folder_symmetric '\activity_*.npy']);
As_symmetric = {act_files.name};
network_files = dir([data_folder_symmetric '\newW_*.npy']);
Ws_symmetric = {network_files.name};
symmetric = get_model_data(Ws_symmetric,As_symmetric,data_folder_symmetric);


%%

figure
subplot(1,2,1);
hold on
scatter(unconstrained.performance*0,unconstrained.performance,'ok','filled')
scatter(unconstrained.performance*0,mean(unconstrained.performance),'dr','filled')
scatter(symmetric.performance*0+1,symmetric.performance,'ok','filled')
scatter(symmetric.performance*0+1,mean(symmetric.performance),'dr','filled')
hold off
xlim([-1 2])
ylim([0 350])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'unconstrained', 'symmetric', ''};
ylabel('performance (sum of absolute difference between model and data)')

subplot(1,2,2);
hold on
scatter(unconstrained.CminusR*0,unconstrained.CminusR,'ok','filled')
scatter(unconstrained.CminusR*0,mean(unconstrained.CminusR),'dr','filled')
scatter(symmetric.CminusR*0+1,symmetric.CminusR,'ok','filled')
scatter(symmetric.CminusR*0+1,mean(symmetric.CminusR),'dr','filled')
yline(0)
hold off
xlim([-1 2])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'unconstrained', 'symmetric', ''};
ylabel('CtoR - RtoC (negative means more RFA to CFA connection strength')


%% idk what we should do to measure the activity in the different regions/regimes

u_from_c = squeeze(mean(unconstrained.no_stim_from_c,2));
u_from_r = squeeze(mean(unconstrained.no_stim_from_r,2));

s_from_c = squeeze(mean(symmetric.no_stim_from_c,2));
s_from_r = squeeze(mean(symmetric.no_stim_from_r,2));

% now just get the mean activity from 0 to 50ms after the reach onset
window = 25:75;
mean_act_UfC = mean(u_from_c(window,:),1);
mean_act_UfR = mean(u_from_r(window,:),1);

mean_act_SfC = mean(s_from_c(window,:),1);
mean_act_SfR = mean(s_from_r(window,:),1);

%%

str_window = [window(1) window(end)]-50;
sz = 10;
figure;
subplot(1,2,1)
hold on
scatter(mean_act_UfC*0,mean_act_UfC,'ob','filled')
scatter(0,mean(mean_act_UfC),sz*10,'dk','filled')
scatter(mean_act_UfR*0+1,mean_act_UfR,'or','filled')
scatter(1,mean(mean_act_UfR),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([-0.3 0.4])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'CFA neurons projecting to RFA', 'RFA neurons projecting to CFA', ''};
ylabel(['mean activity of inter regional projecting neurons in window ' num2str(str_window) ' after reach onset'] )
title('unconstrained (each dot is a model)')

subplot(1,2,2)
hold on
scatter(mean_act_SfC*0,mean_act_SfC,'ob','filled')
scatter(0,mean(mean_act_SfC),sz*10,'dk','filled')
scatter(mean_act_SfR*0+1,mean_act_SfR,'or','filled')
scatter(1,mean(mean_act_SfR),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([-0.3 0.4])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'CFA neurons projecting to RFA', 'RFA neurons projecting to CFA', ''};
ylabel(['mean activity of inter regional projecting neurons in window ' num2str(str_window) ' after reach onset'] )
title('symmetric (each dot is a model)')



%%
% figure;
% subplot(2,1,1);
% hold on
% h1 = plot(u_from_c,'-b');
% h2 = plot(u_from_r,'-r');
% hold off
% title('unconstrained mean activity of inter regional neurons')
% legend([h1(1) h2(2)],{'CFA neurons projecting to RFA','RFA neurons projecting to CFA'})
% 
% subplot(2,1,2);
% hold on
% h1 = plot(s_from_c,'-b');
% h2 = plot(s_from_r,'-r');
% hold off
% title('symmetric mean activity of inter regional neurons')
% legend([h1(1) h2(2)],{'CFA neurons projecting to RFA','RFA neurons projecting to CFA'})
% 
% sgtitle('each line is one model')

%% intra regional stuff
sz = 15;

figure;
subplot(1,2,1);
hold on
scatter(unconstrained.avg_CtoC_exc,unconstrained.var_CtoC_exc,sz,'ob','filled')
scatter(unconstrained.avg_RtoR_exc,unconstrained.var_RtoR_exc,sz,'or','filled')

scatter(symmetric.avg_CtoC_exc,symmetric.var_CtoC_exc,sz,'+b')
scatter(symmetric.avg_RtoR_exc,symmetric.var_RtoR_exc,sz,'+r')
hold off

xlabel('average connection strength')
ylabel('variance of connection strengths')
legend('CFA unconstrained','RFA unconstrained','CFA symmetric','RFA symmetric')
title('excitatory')

subplot(1,2,2);
hold on
scatter(unconstrained.avg_CtoC_inh,unconstrained.var_CtoC_inh,sz,'ob','filled')
scatter(unconstrained.avg_RtoR_inh,unconstrained.var_RtoR_inh,sz,'or','filled')

scatter(symmetric.avg_CtoC_inh,symmetric.var_CtoC_inh,sz,'+b')
scatter(symmetric.avg_RtoR_inh,symmetric.var_RtoR_inh,sz,'+r')
hold off

xlabel('average connection strength')
ylabel('variance of connection strengths')
legend('CFA unconstrained','RFA unconstrained','CFA symmetric','RFA symmetric')
title('inhibitory')

%% same as above but grouping differently
figure;
subplot(1,2,1);
hold on
scatter(unconstrained.avg_CtoC_exc,unconstrained.avg_CtoC_inh,sz,'ob','filled')
scatter(unconstrained.avg_RtoR_exc,unconstrained.avg_RtoR_inh,sz,'or','filled')
scatter(symmetric.avg_CtoC_exc,symmetric.avg_CtoC_inh,sz,'+b')
scatter(symmetric.avg_RtoR_exc,symmetric.avg_RtoR_inh,sz,'+r')
refline(-2.3,0.03801)
hold off
xlabel('average exc connection strength')
ylabel('average inh connection strength')
legend('CFA unconstrained','RFA unconstrained','CFA symmetric','RFA symmetric')
title('averages')

subplot(1,2,2);
hold on
scatter(unconstrained.var_CtoC_exc,unconstrained.var_CtoC_inh,sz,'ob','filled')
scatter(unconstrained.var_RtoR_exc,unconstrained.var_RtoR_inh,sz,'or','filled')
scatter(symmetric.var_CtoC_exc,symmetric.var_CtoC_inh,sz,'+b')
scatter(symmetric.var_RtoR_exc,symmetric.var_RtoR_inh,sz,'+r')
hold off
xlabel('variance of exc connection strength')
ylabel('variance of inh connection strength')
legend('CFA unconstrained','RFA unconstrained','CFA symmetric','RFA symmetric')
title('variances')



%%
sz = 5;
figure;
subplot(2,2,1)
hold on
scatter(unconstrained.avg_CtoC_exc*0,unconstrained.avg_CtoC_exc,sz,'ob','filled')
scatter(0,mean(unconstrained.avg_CtoC_exc),sz*10,'dk','filled')
scatter(unconstrained.avg_RtoR_exc*0+1,unconstrained.avg_RtoR_exc,sz,'or','filled')
scatter(1,mean(unconstrained.avg_RtoR_exc),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([0.06 0.095])
ax = gca;
ax.XTick = [-1,0,1.2];
ax.XTickLabels = {'', 'CFA', 'RFA',''};
ylabel('Average excitatory synapse strength')
title('unconstrained (each dot is a model)')
legend('CFA to CFA','','RFA to RFA')

subplot(2,2,2);
hold on
scatter(symmetric.avg_CtoC_exc*0,symmetric.avg_CtoC_exc,sz,'ob','filled')
scatter(0,mean(symmetric.avg_CtoC_exc),sz*10,'dk','filled')
scatter(symmetric.avg_RtoR_exc*0+1,symmetric.avg_RtoR_exc,sz,'or','filled')
scatter(1,mean(symmetric.avg_RtoR_exc),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([0.06 0.095])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'CFA', 'RFA',''};
ylabel('Average excitatory synapse strength')
title('symmetric (each dot is a model)')

subplot(2,2,3);
hold on
scatter(unconstrained.avg_CtoC_inh*0+0,unconstrained.avg_CtoC_inh,sz,'ob','filled')
scatter(0,mean(unconstrained.avg_CtoC_inh),sz*10,'dk','filled')
scatter(unconstrained.avg_RtoR_inh*0+1,unconstrained.avg_RtoR_inh,sz,'or','filled')
scatter(1,mean(unconstrained.avg_RtoR_inh),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([-0.2 -0.1])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'CFA', 'RFA',''};
ylabel('Average inhibitory synapse strength')

subplot(2,2,4);
hold on
scatter(symmetric.avg_CtoC_inh*0,symmetric.avg_CtoC_inh,sz,'ob','filled')
scatter(0,mean(symmetric.avg_CtoC_inh),sz*10,'dk','filled')
scatter(symmetric.avg_RtoR_inh*0+1,symmetric.avg_RtoR_inh,sz,'or','filled')
scatter(1,mean(symmetric.avg_RtoR_inh),sz*10,'dk','filled')
hold off
xlim([-1 2])
ylim([-0.2 -0.1])
ax = gca;
ax.XTick = [-1,0,1,2];
ax.XTickLabels = {'', 'CFA', 'RFA',''};
ylabel('Average inhibitory synapse strength')   

%%




Q= get_model_data(Ws_unconstrained,As_unconstrained,data_folder_unconstrained);


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
