%% load the connectivity matrix and do stuff with it

addpath(genpath('Z:/code'));
data_folder = uigetdir();

%%
files = dir([data_folder '\newW_*.npy']);
Ws = {files.name};

W = zeros(1000,1000,length(Ws));
CtoR = [];
RtoC = [];
differences = [];
performance = [];
for i=1:length(Ws)
    w_file = Ws{i};
    perf = split(w_file,'_');
    perf = perf{end};
    perf = perf(1:end-4);   % get rid of the .npy
    perf = str2num(perf);
    performance = [performance; perf];
    %     w = readNPY(['outputweights_symmetric/' w_file]);
    w = readNPY([data_folder '\' w_file]);
    W(:,:,i) = w;
    cfa_to_rfa = w(501:end,1:500);
    rfa_to_cfa = w(1:500,501:end);
    d = sum(rfa_to_cfa(:)) - sum(cfa_to_rfa(:));

    CtoR = [CtoR; cfa_to_rfa(:)];
    RtoC = [RtoC; rfa_to_cfa(:)];
    differences = [differences; d];
end
CtoR(CtoR==0) = [];
RtoC(RtoC==0) = [];

% sum(RtoC)
% sum(CtoR)
%
figure;
scatter(performance,differences,'.k')
xlabel('performance')
ylabel('difference RtoC - CtoR')

%%
% RtoC_max = max(RtoC);
% CtoR_max = max(CtoR);
% lowest_max = min(RtoC_max,CtoR_max);
% 
% RtoC(RtoC>lowest_max) = nan;
% CtoR(CtoR>lowest_max) = nan;


%%
% two sample kstest
[h,p] = kstest2(CtoR,RtoC)


bw = 0.025;
figure;
subplot(3,1,1);
hold on
histogram(CtoR,'BinWidth',bw);
histogram(RtoC,'BinWidth',bw);
hold off
legend({'cfa->rfa','rfa->cfa'})

subplot(3,1,2)
hold on
boxplot([RtoC CtoR],'orientation', 'horizontal', 'Whisker',Inf);
hold off
annotation('textbox',[.2 .27 .3 .3],'String','cfa->rfa','FitBoxToText','on');
annotation('textbox',[.2 .17 .3 .3],'String','rfa->cfa','FitBoxToText','on');

subplot(3,1,3);
hold on
cdfplot(CtoR);
cdfplot(RtoC);
hold off
legend({'cfa->rfa','rfa->cfa'})
title(sprintf('2 sample ks test p value = %0.11f',p));


%%










%%