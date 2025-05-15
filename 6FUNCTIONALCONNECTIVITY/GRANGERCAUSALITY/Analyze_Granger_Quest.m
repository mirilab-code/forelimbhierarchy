% ... and now let's plot the results
date = '03302020';
str1 = sprintf('C:/Users/mirilab/Documents/Adam/reaching/10k_granger_trains/%s/%s_Out.mat',...
    date,date);
plotResults_Granger_new(str1);

%%
load(str1)
psi = OutStruct.Psi2;
num_CFA = size(psi,1)/2;
conns = zeros(1,4);
conns(1) = length(find(psi(1:num_CFA,1:num_CFA)));
conns(2) = length(find(psi(1:num_CFA,num_CFA+1:end)));
conns(3) = length(find(psi(num_CFA+1:end,1:num_CFA)));
conns(4) = length(find(psi(num_CFA+1:end,num_CFA+1:end)));
all_str = ['CFA';'CFA';'RFA';'CFA';'CFA';'RFA';'RFA';'RFA'];
str = sprintf('\n');
disp(str)
disp('Self_Connections Included')
for ii = 1:4
    str = sprintf('%s->%s: %d',all_str(2*ii-1,:),...
            all_str(2*ii,:),conns(ii));
    disp(str)
end
conns_noself = conns;
for ii = 1:size(psi,1)
    if psi(ii,ii) == 1 && ii <= num_CFA
        conns_noself(1) = conns_noself(1) - 1;
    elseif psi(ii,ii) == 1 && ii > num_CFA
        conns_noself(4) = conns_noself(4) - 1;
    end
end
str = sprintf('\n');
disp(str)
disp('Self_Connections Removed')
for ii = 1:4
    str = sprintf('%s->%s: %d',all_str(2*ii-1,:),...
            all_str(2*ii,:),conns_noself(ii));
    disp(str)
end
%%
D = OutStruct.glm_dev_ratio;
ht = OutStruct.neuronHistoryRegressorNBins;
nNeurons = size(psi,1);
P = zeros(nNeurons,nNeurons);
newpsi = zeros(nNeurons,nNeurons);
for n = 1:nNeurons
     P(n,:) = 1 - chi2cdf(D(n,:),ht(n));
end
CFA_CFAp = P(1:num_CFA,1:num_CFA);
CFA_RFAp = P(num_CFA+1:end,1:num_CFA);
RFA_CFAp = P(1:num_CFA,num_CFA+1:end);
RFA_RFAp = P(num_CFA+1:end,num_CFA+1:end);
sep_p = {CFA_CFAp,CFA_RFAp,RFA_CFAp,RFA_RFAp};
figure;
histogram(P(:),40)
xlabel('p-value')
str = sprintf('%s ALL p-values',date');
title(str)
figure;
subplot(2,2,1)
for ii = 1:4
    subplot(2,2,ii);
    histogram(sep_p{ii}(:),40);
    str = sprintf('%s->%s',all_str(2*ii-1,:),all_str(2*ii,:));
    title(str);
end

% figure;
% imagesc(P)
% figure;
% histogram(D(:))
figure;
imagesc(D)
str = sprintf('%s: Difference of Model Deviations',date);
colorbar
title(str)
figure;
scatter(1:40,ht)
xlabel('Neuron ID')
ylabel('Number of History Bins')
xline(20.5);
str = sprintf('%s: History Bins for Each Neuron',date);
title(str);

% for ii = 1:nNeurons*nNeurons
%    if P(ii)<0.05
%        newpsi(ii) = 1;
%    end
% end
% figure;
% imagesc(newpsi)
% figure;
% imagesc(OutStruct.Psi1)

%%
histBins = OutStruct.historyRegressorNBins;
neurBins = OutStruct.neuronHistoryRegressorNBins;
D_org = cell(1,5);
for ii=1:length(histBins)
    temp = find(neurBins==histBins(ii));
    D_org{ii} = D(temp,:);
    D_org{ii} = D_org{ii}(:);
    figure;
    edges = 0:2:50;
    histogram(D_org{ii},edges);
    str = sprintf('Number of History Bins = %d',histBins(ii));
    title(str)
    x = 1:0.1:50;
    y = chi2pdf(x,histBins(ii));
    hold on
    plot(x,length(D_org{ii})*y,'r');
end

binMeansAdj = zeros(1,5);
for ii = 1:length(histBins)
    binMeansAdj(ii) = mean(D_org{ii})/histBins(ii);
end
disp(binMeansAdj)

%%
