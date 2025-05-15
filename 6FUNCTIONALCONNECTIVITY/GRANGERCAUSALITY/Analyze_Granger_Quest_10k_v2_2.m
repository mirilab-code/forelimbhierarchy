% ... and now let's plot the results

dates = ['03272020';'03302020';'03312020';'04012020';'04022020';'04032020';...
    '09222020';'09232020';'09302020';'10012020';'10022020';'10052020';...
    '10062020';'10242020';'10262020';'10272020';'10282020';'11102020';...
    '11122020';'11152020';'11172020'];
dates_run = ['03272020';'03302020';'03312020';'04012020';'04022020';'04032020';...
    '09222020';'09232020';'09302020';'10012020';'10022020';'10052020';...
    '10062020';'10242020';'10262020';'10272020';'10282020';'11102020';...
    '11122020';'11152020';'11172020'];

all_pvals_C2C = [];
all_pvals_R2R = [];
all_pvals_R2C = [];
all_pvals_C2R = [];
all_D_C2C = [];
all_D_R2R = [];
all_D_R2C = [];
all_D_C2R = [];

for session_num = 1:size(dates_run,1)
    date = dates_run(session_num,:);
    ii = 1;
    while ii < length(dates(:,1))+1
        if strcmp(dates(ii,:),date)
            date_number = ii;
            break
        end
        ii = ii + 1;
    end

    str = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_num_neurons.mat',...
        date,date);
    load(str);
    str = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_final_list.mat',...
        date,date);
    load(str);
    str = sprintf('Z://10k_FRpairs//top_ten_thou_upd.mat');
    load(str);
    num_CFA = num_neurons(1);
    num_RFA = num_neurons(2);
    total_neurons = num_CFA + num_RFA;
    session_inds = find(top_ten_thou_upd(:,1)==date_number);
    CFA_RFA_pairs = top_ten_thou_upd(session_inds,3:4);
    maxCFApairInd = min(find(final_list(:,1)>num_CFA))-1;
    maxRFApairInd = 2*maxCFApairInd;
    maxR2CpairInd = maxRFApairInd+length(session_inds);
    CFA_unique_inds = sort(unique(CFA_RFA_pairs(:,1)));
    RFA_unique_inds = sort(unique(CFA_RFA_pairs(:,2)));

    if ismember(date_number,[1,15])
        str1 = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_Out_fixed_only_parneur.mat',...
            date,date);
        load(str1);
    elseif date_number == 18
        str1 = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_Out_fixed_only_parneur_bottomhalf.mat',...
            date,date);
        load(str1)
    else
        str1 = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_Out_fixed_only.mat',...
            date,date);
        load(str1);
    end

%     str1 = sprintf('C:/Users/mirilab/OneDrive - Northwestern University/Documents/Adam/reaching/10k_granger_trains/%s/%s_Out_fixed_only_parneur_pre.mat',...
%         date,date);
%     load(str1);
 
    Psi1_old = OutStruct.Psi1;
    P = zeros(nNeurons,nNeurons);
    D = OutStruct.glm_dev_ratio;
    ht = OutStruct.neuronHistoryRegressorNBins;
    for ii = 1:nNeurons
        for jj= 1:nNeurons
            if ~ismember([ii,jj],final_list,'rows')
                D(ii,jj) = NaN;
            end
            if ismember(ii,neural_warnings{date_number})
                D(ii,:) = NaN;
            end
        end
    end
    if ~isempty(neural_pair_warnings{date_number})
        for ii = 1:size(neural_pair_warnings{date_number},1)
            temp1 = neural_pair_warnings{date_number}(ii,1);
            temp2 = neural_pair_warnings{date_number}(ii,2);
            D(temp1,temp2) = NaN;
        end
    end

    figure;
    imagesc(Psi1_old)
    str = sprintf('%s: Recovered Connectivity (Before FDR)',date);
    title(str)
    xlabel('Source')
    ylabel('Target')
    conns = zeros(1,4);
    conns(1) = length(find(Psi1_old(1:num_CFA,1:num_CFA)));
    conns(2) = length(find(Psi1_old(1:num_CFA,num_CFA+1:end)));
    conns(3) = length(find(Psi1_old(num_CFA+1:end,1:num_CFA)));
    conns(4) = length(find(Psi1_old(num_CFA+1:end,num_CFA+1:end)));
    all_str = ['CFA';'CFA';'RFA';'CFA';'CFA';'RFA';'RFA';'RFA'];
    fprintf('\n');
    disp(date)
    for ii = 1:4
        fprintf('%s->%s: %d\n',all_str(2*ii-1,:),...
                all_str(2*ii,:),conns(ii));
    end
    disp('Intra-Pairs = ')
    disp(maxCFApairInd)
    disp('Inter-Pairs = ')
    disp(maxR2CpairInd-maxRFApairInd)
    figure;
    imagesc(D)
    if strcmp(date,'10262020')
        disp(D)
    end
 
    str = sprintf('%s: Difference in Model Deviation',date);
    colorbar
    title(str)
    xline(num_CFA+0.5,'r');
    yline(num_CFA+0.5,'r');
    xlabel('Source Neuron')
    ylabel('Target Neruon')
    for n = 1:total_neurons
         P(n,:) = 1 - chi2cdf(D(n,:),ht(n));
    end
    figure;
    imagesc(P)
    str = sprintf('%s: p-Values',date);
    colorbar
    title(str)
    xline(num_CFA+0.5,'r');
    yline(num_CFA+0.5,'r');
    xlabel('Source Neuron')
    ylabel('Target Neruon')
    for temp = 1:maxCFApairInd
        all_pvals_C2C = [all_pvals_C2C,P(final_list(temp,1),final_list(temp,2))];
        all_D_C2C = [all_D_C2C,D(final_list(temp,1),final_list(temp,2))];
    end
    for temp = maxCFApairInd+1:maxRFApairInd
        all_pvals_R2R = [all_pvals_R2R,P(final_list(temp,1),final_list(temp,2))];
        all_D_R2R = [all_D_R2R,D(final_list(temp,1),final_list(temp,2))];
    end
    for temp = maxRFApairInd+1:maxR2CpairInd
        all_pvals_R2C = [all_pvals_R2C,P(final_list(temp,1),final_list(temp,2))];
        all_D_R2C = [all_D_R2C,D(final_list(temp,1),final_list(temp,2))];
    end
    for temp = maxR2CpairInd+1:size(final_list(:,1))
        all_pvals_C2R = [all_pvals_C2R,P(final_list(temp,1),final_list(temp,2))];
        all_D_C2R = [all_D_C2R,D(final_list(temp,1),final_list(temp,2))];
    end
    
    pairs_square = zeros(total_neurons);
    pairs_square_C2R = zeros(total_neurons);
    pairs_square_R2C = zeros(total_neurons);
    Psi1_old = OutStruct.Psi1;
    Psi2_old = OutStruct.Psi2;
    Psi3 = zeros(total_neurons);
    for ii = 1:size(CFA_RFA_pairs,1)
        CFA = find(CFA_unique_inds == CFA_RFA_pairs(ii,1));
        RFA = find(RFA_unique_inds == CFA_RFA_pairs(ii,2))+num_CFA;
        pairs_square(CFA,RFA) = 1;
        pairs_square_C2R(CFA,RFA) = 1;
        pairs_square(RFA,CFA) = 1;
        pairs_square_R2C(RFA,CFA) = 1;
        if Psi1_old(CFA,RFA)==1
            Psi3(CFA,RFA)=1;
        end
        if Psi1_old(RFA,CFA)==1
            Psi3(RFA,CFA)=1;
        end
    end
end
%%
all_D_both = [all_D_C2C,all_D_R2R,all_D_R2C,all_D_C2R];
all_D_C2R(isnan(all_D_C2R)) = [];
all_D_R2C(isnan(all_D_R2C)) = [];
all_D_C2C(isnan(all_D_C2C)) = [];
all_D_R2R(isnan(all_D_R2R)) = [];

pairs_diff = length(all_D_R2C)-length(all_D_C2R);
if pairs_diff>0
    extra_pairs_inds = randperm(length(all_D_R2C),pairs_diff);
    all_D_R2C(extra_pairs_inds) = [];
elseif pairs_diff<0
    extra_pairs_inds = randperm(length(all_D_C2R),abs(pairs_diff));
    all_D_C2R(extra_pairs_inds) = [];
end


figure;
edges = 0:1:(10*ceil(max(all_D_both)/10));
histogram(all_D_C2R,edges)
hold on
histogram(all_D_R2C,edges)
histogram(all_D_C2C,edges)
histogram(all_D_R2R,edges)
xlabel('Difference of Model Deviations')
ylabel('Number of Pairs')
legend('CFA -> RFA','RFA -> CFA');
title('Connectivity Strength - 10k Pairs Only')
C2Rmed = median(all_D_C2R);
R2Cmed = median(all_D_R2C);
% xline(C2Rmed,'b')
% xline(R2Cmed,'r')

all_pvals_C2R(isnan(all_pvals_C2R)) = [];
all_pvals_R2C(isnan(all_pvals_R2C)) = [];

if pairs_diff>0
    all_pvals_R2C(extra_pairs_inds) = [];
elseif pairs_diff<0
    all_pvals_C2R(extra_pairs_inds) = [];
end
figure;
edges = 0:0.025:1;
histogram(all_pvals_C2R,edges)
hold on
histogram(all_pvals_R2C,edges)
% histogram(all_pvals_C2C,edges)
% histogram(all_pvals_R2R,edges)
xlabel('p-value')
ylabel('Number of Pairs')
legend('CFA -> RFA','RFA -> CFA');
title('P-Values of Difference of Model Deviation - 10k Pairs Only')

C2R_signif = find(all_pvals_C2R<0.05);
R2C_signif = find(all_pvals_R2C<0.05);

C2R_D_signif = all_D_C2R(C2R_signif);
R2C_D_signif = all_D_R2C(R2C_signif);

figure;
edges = 18:1:(10*ceil(max(all_D_both)/10));
histogram(C2R_D_signif,edges)
%set(gca,'xscale','log')
set(gca,'yscale','log')
hold on
histogram(R2C_D_signif,edges)
xlabel('Difference of Model Deviations')
ylabel('Number of Pairs')
legend('CFA -> RFA','RFA -> CFA');
title('Connectivity Strength - 10k Pairs Only')

%%
[fdr_c2r,q_c2r,priori_c2r] = mafdr(all_pvals_C2R,'Method','bootstrap','Showplot',true);
[fdr_r2c,q_r2c,priori_r2c] = mafdr(all_pvals_R2C,'Method','bootstrap','Showplot',true);
b_c2r= find (fdr_c2r<.2);
b_r2c= find (fdr_r2c<.2);
C2R_sig = all_D_C2R(b_c2r);
R2C_sig = all_D_R2C(b_r2c);

figure;
histogram(C2R_sig,'BinWidth',5)
hold on
histogram(R2C_sig,'BinWidth',5)
legend('CFA -> RFA','RFA -> CFA');
xlabel('Difference of Model Deviations')
ylabel('Number of Pairs')
title('Connectivity Strength (FDR < 0.2)')

%%
%%%

pairs_diff = length(all_D_C2C)-length(all_D_R2R);
if pairs_diff>0
    extra_pairs_inds = randperm(length(all_D_C2C),pairs_diff);
    all_D_C2C(extra_pairs_inds) = [];
elseif pairs_diff<0
    extra_pairs_inds = randperm(length(all_D_R2R),abs(pairs_diff));
    all_D_R2R(extra_pairs_inds) = [];
end

figure;
edges = 0:1:(10*ceil(max(all_D_both)/10));
%histogram(all_D_C2R,edges)
hold on
%histogram(all_D_R2C,edges)
histogram(all_D_C2C,edges)
histogram(all_D_R2R,edges)
xlabel('Difference of Model Deviations')
ylabel('Number of Pairs')
legend('CFA -> CFA','RFA -> RFA');
title('Connectivity Strength - Intra Pairs Only')
C2Cmed = median(all_D_C2C);
R2Rmed = median(all_D_R2R);
% xline(C2Rmed,'b')
% xline(R2Cmed,'r')

all_pvals_C2C(isnan(all_pvals_C2C)) = [];
all_pvals_R2R(isnan(all_pvals_R2R)) = [];
if pairs_diff>0
    all_pvals_C2C(extra_pairs_inds) = [];
elseif pairs_diff<0
    all_pvals_R2R(extra_pairs_inds) = [];
end

figure;
edges = 0:0.025:1;
%histogram(all_pvals_C2R,edges)
hold on
%histogram(all_pvals_R2C,edges)
histogram(all_pvals_C2C,edges)
histogram(all_pvals_R2R,edges)
xlabel('p-value')
ylabel('Number of Pairs')
legend('CFA -> CFA','RFA -> RFA');
title('P-Values of Difference of Model Deviation - intra Pairs Only')

% C2C_signif = find(all_pvals_C2C<0.05);
% R2C_signif = find(all_pvals_R2C<0.05);
% 
% C2R_D_signif = all_D_C2R(C2R_signif);
% R2C_D_signif = all_D_R2C(R2C_signif);
% 
% figure;
% edges = 18:1:(10*ceil(max(all_D_both)/10));
% histogram(C2R_D_signif,edges)
% %set(gca,'xscale','log')
% set(gca,'yscale','log')
% hold on
% histogram(R2C_D_signif,edges)
% xlabel('Difference of Model Deviations')
% ylabel('Number of Pairs')
% legend('CFA -> RFA','RFA -> CFA');
% title('Connectivity Strength - 10k Pairs Only')

%%
figure;
[fdr_c2c,q_c2c,priori_c2c] = mafdr(all_pvals_C2C,'Method','bootstrap','Showplot',true);
[fdr_r2r,q_r2r,priori_r2r] = mafdr(all_pvals_R2R,'Method','bootstrap','Showplot',true);
b_c2c= find (fdr_c2c<.2);
b_r2r= find (fdr_r2r<.2);
C2C_sig = all_D_C2C(b_c2c);
R2R_sig = all_D_R2R(b_r2r);

figure;
histogram(C2C_sig,'BinWidth',5)
hold on
histogram(R2R_sig,'BinWidth',5)
legend('CFA -> CFA','RFA -> RFA');
xlabel('Difference of Model Deviations')
ylabel('Number of Pairs')
title('Connectivity Strength (FDR < 0.2)')