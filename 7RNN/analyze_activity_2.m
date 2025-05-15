

%%
mean_sym_act_from_c = squeeze(mean(symmetric.no_stim_from_c,2));
mean_sym_act_from_r = squeeze(mean(symmetric.no_stim_from_r,2));

mean_unc_act_from_c = squeeze(mean(unconstrained.no_stim_from_c,2));
mean_unc_act_from_r = squeeze(mean(unconstrained.no_stim_from_r,2));






figure;
subplot(2,1,1);
hold on
hc = plot(mean_unc_act_from_c,'Color',[0 0 1 0.7]);
plot(mean(mean_unc_act_from_c,2),'-k','LineWidth',2);
hr = plot(mean_unc_act_from_r,'Color',[1 0 0 0.7]);
plot(mean(mean_unc_act_from_r,2),'-k','LineWidth',2);
title('mean unconstrained')


subplot(2,1,2);
hold on
hc = plot(mean_sym_act_from_c,'Color',[0 0 1 0.7]);
plot(mean(mean_sym_act_from_c,2),'-k','LineWidth',2);
hr = plot(mean_sym_act_from_r,'Color',[1 0 0 0.7]);
plot(mean(mean_sym_act_from_r,2),'-k','LineWidth',2);
legend([hc(1) hr(1)],{'from CFA == to RFA','from RFA == to CFA'})
title('mean symmetric')


%%
figure;
tiledlayout(5,6);
for i=1:30
    nexttile;
    hold on

%     quants = [0.75 0.5 0.25];
%     q_from_c = quantile(from_c(:,:,i)',quants);
%     q_from_r = quantile(from_r(:,:,i)',quants);
%     hold on;
%     hc = plot(q_from_c','Color',[0 0 1 0.9]);
%     hr = plot(q_from_r','Color',[1 0 0 0.9]);

    a = mean_sym_act_from_c(:,i);
    b = mean_sym_act_from_r(:,i);
    plot(a,'-b')
    plot(b,'-r')
%     plot(a+b,'-k')
%     plot(mean([a b],2),'-k');

end
legend([hc(1) hr(1)],{'from CFA','from RFA'})
sgtitle('symmetric: 25 50 75 quantiles of RFA/CFA activity in projecting neurons weighted by their connection strength')


%%
plot(mean_sym_act_from_c + mean_sym_act_from_r,'-k')
title('from cfa + from rfa')


%%
quants_s_from_c = quantile(symmetric.no_stim_from_c,[0.75 0.5 0.25],2);
quants_s_from_r = quantile(symmetric.no_stim_from_r,[0.75 0.5 0.25],2);


%%
% subplot(5,1,3);
% hold on
% plot(squeeze(quants_s_from_c(:,1,:)),'Color',[0 0 1 0.7])
% plot(squeeze(quants_s_from_r(:,1,:)),'Color',[1 0 0 0.7])
% title('nth quantile')
% 
% subplot(5,1,4);
% hold on
% plot(squeeze(quants_s_from_c(:,3,:)),'Color',[0 0 1 0.7])
% plot(squeeze(quants_s_from_r(:,3,:)),'Color',[1 0 0 0.7])
% title('nth quantile')
%%

modn = modn+1;

quants = [0.75 0.5 0.25];
q_from_c = quantile(from_c(:,:,modn)',quants);
q_from_r = quantile(from_r(:,:,modn)',quants);

%


hold on;
hc = plot(q_from_c','Color',[0 0 1 0.9]);
hr = plot(q_from_r','Color',[1 0 0 0.9]);
legend([hc(1) hr(1)],{'from CFA','from RFA'})

%%
% mm = 2;
mm = mm+1;
figure(1)
clf("reset")
hold on
plot(from_c(:,:,mm),'Color',[0 0 1 0.2])
plot(from_r(:,:,mm),'Color',[1 0 0 0.2])