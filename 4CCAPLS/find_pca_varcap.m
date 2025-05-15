%%

target_vc = 0.95;

totalvarcap_RFA = {};
totalvarcap_CFA = {};
for i=1:length(D)
% for i=6:6
    disp(i);
    A = D(i);
    [~,~,latentRFA] = pca(A.RFA_RG');
    [~,~,latentCFA] = pca(A.CFA_RG(:,:,lag0)');

    latentRFA = latentRFA/sum(latentRFA);
    latentCFA = latentCFA/sum(latentCFA);
    varcapRFA = cumsum(latentRFA);
    varcapCFA = cumsum(latentCFA);

    totalvarcap_RFA{i} = varcapRFA;
    totalvarcap_CFA{i} = varcapCFA;



end
%%
mindimRFA = min(cellfun(@length, totalvarcap_RFA));
mindimCFA = min(cellfun(@length, totalvarcap_CFA));

vc_RFA = [];
vc_CFA = [];

for i=1:length(totalvarcap_CFA)
    vc_rfa = totalvarcap_RFA{i}(1:mindimRFA)';
    vc_RFA = [vc_RFA; vc_rfa];
    vc_cfa = totalvarcap_CFA{i}(1:mindimCFA)';
    vc_CFA = [vc_CFA; vc_cfa];
end

%%
mean_vc_RFA = mean(vc_RFA);
mean_vc_CFA = mean(vc_CFA);

nRFA = find(mean_vc_RFA>target_vc,1)
nCFA = find(mean_vc_CFA>target_vc,1)























%%