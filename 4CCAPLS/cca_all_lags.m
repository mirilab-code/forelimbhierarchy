function R = cca_all_lags(M1,M2_lags,lags,nPCs)

R = {};
vc_RFA = {};
vc_CFA = {};


for i=1:length(lags)
    lag = lags(i);
    [~,~,r,~,~,varcap1,varcap2] = pca_then_cca(M1,M2_lags(:,:,i),nPCs);
    R{i} = r;

%     figure;
%     hold on
%     plot(r, 'Color',[0 0 0 0.3])

    vc_RFA{i} = varcap1.cca;
    vc_CFA{i} = varcap2.cca;
end


weighted_Rs = [];
Rs = [];
first_comp_R = [];


for i=1:length(R)
    w = mean([vc_RFA{i} vc_CFA{i}],2);
    wr = dot(w,R{i});
    first_comp_R = [first_comp_R; R{i}(1)];
    r = mean(R{i});
    weighted_Rs = [weighted_Rs; wr];
    Rs = [Rs; r];


end

R = struct;
R.weighted = weighted_Rs';
R.unweighted = Rs';
R.first_comp = first_comp_R';
R.varcapRFA = vc_RFA;
R.varcapCFA = vc_CFA;

end