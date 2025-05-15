%% just the PCs



ncomps = num_PCs;

PCA_varcapRFA = [];
PCA_varcapCFA = [];
CCA_varcapRFA = [];
CCA_varcapCFA = [];

CORRS = [];
for i=1:1
    disp(i);
    lags = D(i).lags;
    lag0 = find(lags==0);
    RFA = D(i).RFA_RG;
    CFA = D(i).CFA_RG(:,:,lag0);
    

    [coefRFA,scoreRFA,latentRFA] = pca(RFA');
    [coefCFA,scoreCFA,latentCFA] = pca(CFA');
    vRFA = latentRFA / sum(latentRFA);
    vCFA = latentCFA / sum(latentCFA);

    vRFA = vRFA(1:ncomps);
    vCFA = vCFA(1:ncomps);

end

%%
x = -99:100;

figure;
for i=1:2
    for j=1:8
        ind = sub2ind([8 2], j, i);

        subplot(4,8,ind)
        window = 1+(j-1)*200:((j-1)*200)+200;
        disp([window(1) window(end)])
        
        trace = scoreRFA(window,i);
        plot(x,trace)
        ylim([-0.15 0.15])
%         pause(0.5);
        if(j==1)
            vc_rfa = vRFA(i);
            ylabel(sprintf('RFA PC %d,\n varcap: %0.3f ',i,vc_rfa))
        end
    end
end
stop = ind;

for i=1:2
    for j=1:8
        ind = sub2ind([8 2], j, i);

        subplot(4,8,stop+ind)
        window = 1+(j-1)*200:((j-1)*200)+200;
        disp([window(1) window(end)])
        
        trace = scoreCFA(window,i);
        plot(x,trace)
        ylim([-0.15 0.25])
%         pause(0.5);
        if(j==1)
            vc_cfa = vCFA(i);
            ylabel(sprintf('CFA PC %d,\n varcap: %0.3f ',i,vc_cfa))
        end
    end
end




















%%