addpath(genpath('Z:\Scripts\boundedline\'));

as = {'a048','a050','a051','ss2','MA1','MA2'};

N = length(as);


data = {};
lag_unweightedRs = {};
lag_weightedRs = {};
for i=1:length(as)
    fname = ['results/' as{i} '_CCA.mat']
    data{i} = load(fname);
    lags = data{1}.lags;
    lag_unweightedRs{i} = data{i}.unweightedRs;
    lag_weightedRs{i} = data{i}.weightedRs;
end

Rs_unweighted = cell2mat(lag_unweightedRs)';
Rs_weighted = cell2mat(lag_weightedRs)';

%% get the mean lag correlations for all the animals.

mean_unweighted = mean(Rs_unweighted);
mean_weighted = mean(Rs_weighted);
SEM_unweighted = std(Rs_unweighted)/sqrt(N);
SEM_weighted = std(Rs_weighted)/sqrt(N);


figure;
hold on
[a,b] = boundedline(lags,mean_unweighted,SEM_unweighted,'-r','alpha');
xlabel('CFA lag (ms)');
ylabel('average correlation');

[a,b] = boundedline(lags,mean_weighted,SEM_weighted,'-b','alpha');

ylim([0.5 1.0]);
legend({'nothing','unweighted','nothing','weighted'});


%% get the correlation for each CV stuff
R_CVs = {};
for i=1:length(as)
    
    thisRs = data{i}.R;
    R0 = thisRs{31};
    R_CVs{i} = R0;
%     minD = min(cellfun(@length,thisRs));
%     thisRs = cellfun(@(x) x(1:minD), thisRs, 'un', 0);
%     R_CVs{i} = cell2mat(thisRs');

end
minD = min(cellfun(@length, R_CVs));
R_CVs = cell2mat(cellfun(@(x) x(1:minD), R_CVs, 'un', 0)');

total_mean_RCVs = mean(R_CVs,1);
total_SEM_RCVs = std(R_CVs)/sqrt(N);


vc = load('results/varcap_total.mat');
vcCFA = vc.varcap_total_CFA;
vcRFA = vc.varcap_total_RFA;
minD = min(min(cellfun(@length, vcCFA)),min(cellfun(@length, vcRFA)))

vcCFA_total = cell2mat(cellfun(@(x) x(1:minD), vcCFA, 'un', 0));
vcRFA_total = cell2mat(cellfun(@(x) x(1:minD), vcRFA, 'un', 0));
mean_varcap_CFA = cumsum(mean(vcCFA_total,2));
SEM_varcap_CFA = std(vcCFA_total,0,2)/sqrt(N);
mean_varcap_RFA = cumsum(mean(vcRFA_total,2));
SEM_varcap_RFA = std(vcRFA_total,0,2)/sqrt(N);
%%
figure;
corr = errorbar(1:minD,total_mean_RCVs,total_SEM_RCVs,'ok');
set(gca,'YColor','k')
ylabel('correlation ');
xlabel('canonical variable');

figure;
hold on
cfa_varcap = errorbar([1:minD]-0.1,mean_varcap_CFA,SEM_varcap_CFA,'og');
rfa_varcap = errorbar([1:minD]+0.1,mean_varcap_CFA,SEM_varcap_RFA,'or');
hold off
legend({'CFA','RFA'})

hold off
xlabel('canonical variable');
ylabel('cumulative percent of variance captured');
ylim([0 1]);
set(gca,'YColor','m')


xlim([0 minD+1]);
legend_data = [corr(1) cfa_varcap(1) rfa_varcap(1)];
legend(legend_data,'correlation','CFA','RFA');
% legend({'nothing','correlation','nothing','CFA','nothing','RFA'});





%% we use the average of these to make the above one^
% figure;
% tiledlayout(2,1);
% 
% nexttile;
% hold on
% for i=1:length(lag_Rs)
%     lags = data{i}.lags;
%     r = lag_unweightedRs{i};
%     [vals,locs] = findpeaks(r);
%     locs = locs - 31;
%     indMax = find(vals==max(vals));
%     locMax = locs(indMax);
%     scatter(locMax,vals(indMax)+0.002, 'vr');
%     
% 
%     plot(lags, r ,'b')
% end
% hold off
% xlabel('CFA lag (ms)');
% ylabel('mean correlation');
% 
% 
% 
% 
% nexttile;
% hold on
% for i=1:length(lag_Rs)
%     lags = data{i}.lags;
%     r = lag_Rs{i};
%     [vals,locs] = findpeaks(r);
%     locs = locs - 31;
%     indMax = find(vals==max(vals));
%     locMax = locs(indMax);
%     scatter(locMax,vals(indMax)+0.002, 'vr');
%     
% 
%     plot(lags, r ,'b')
% end
% hold off
% xlabel('CFA lag (ms)');
% ylabel('mean correlation (weighted by variance capture)');
% 
% 
% 

disp('done!');