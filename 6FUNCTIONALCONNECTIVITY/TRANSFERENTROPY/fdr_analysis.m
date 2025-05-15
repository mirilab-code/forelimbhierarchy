%% false discovery rate

load('Z:\10k_FRpairs\TEoutput\results\TEresults.mat');
%% 

% do the regular FDR
[FDR_CtoR,Q_CtoR,apriori_CtoR] = mafdr(CtoR_P, 'Lambda', 0.5, 'Showplot',false);
[FDR_RtoC,Q_RtoC,apriori_RtoC] = mafdr(RtoC_P, 'Lambda', 0.5, 'Showplot',false);

false_discovery_rate = 0.2;
sig_CtoR = Q_CtoR<false_discovery_rate;
sig_RtoC = Q_RtoC<false_discovery_rate;

nsig_CtoR = sum(sig_CtoR)
nsig_RtoC = sum(sig_RtoC)

%% do the BH FDR
[CtoR_h,~,~,Q_CtoR_bh] = fdr_bh(CtoR_P,false_discovery_rate,'pdep');
[RtoC_h,~,~,Q_RtoC_bh] = fdr_bh(RtoC_P,false_discovery_rate,'pdep');

sig_CtoR_bh = Q_CtoR_bh<false_discovery_rate;
sig_RtoC_bh = Q_RtoC_bh<false_discovery_rate;

nsig_CtoR_bh = sum(sig_CtoR_bh)
nsig_RtoC_bh = sum(sig_RtoC_bh)




%% plotting
figure;
tiledlayout(3,1);
nexttile;
hold on
histogram(RtoC_P, 'BinWidth',0.025);
histogram(CtoR_P, 'BinWidth',0.025);
hold off
legend({'RFA to CFA','CFA to RFA'});
xlabel('p value');
ylabel('count');

nexttile;
hold on
histogram(RtoC_TE(sig_RtoC), 'BinWidth',0.00001);
histogram(CtoR_TE(sig_CtoR), 'BinWidth',0.00001);
hold off
legend({'RFA to CFA','CFA to RFA'});
xlabel('transfer entropy value');
ylabel('count');

nexttile;
B = [CtoR_TE RtoC_TE];
boxplot(B,  'orientation', 'horizontal')
yticklabels({'CFA to RFA','RFA to CFA'});


disp('done!');


%%









