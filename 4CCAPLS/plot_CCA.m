%% plot PLSSVD traces

%% plot traces PLS
lag0 = find(lags==0);
k = 1;
RFAfull = D(k).RFA_RG;
CFAfull = D(k).CFA_RG(:,:,lag0);

dim = 25;

[~,~,r,RFAcca,CFAcca,vcRFA,vcCFA] = pca_then_cca(RFA,CFA,ncomps);


varcapRFA = vcRFA.cca;
varcapCFA = vcCFA.cca;


%%
components = [1 2 5 10];
x = -99:100;

figure;
for i=1:length(components)
    comp = components(i);
    for j=1:8
%         pause(0.25)
        window = 1+(j-1)*200:((j-1)*200)+200;
        disp([window(1) window(end)]);
        ind = sub2ind([8 length(components)], j, i);
        subplot(length(components),8,ind)
        trace1 = RFAcca(window,comp);
        trace2 = CFAcca(window,comp);
    
        hold on
        plot(x,trace1)
        if(j==1)
%             covar = covs(comp);
            this_r = r(comp);
            vc_rfa = varcapRFA(comp);
            vc_cfa = varcapCFA(comp);
            ylabel(sprintf('PLS component %d, \n  covariance: %0.3f \n vcRFA: %0.3f vcCFA: %0.3f ',comp,this_r,vc_rfa,vc_cfa))
        end
%         yyaxis right
        plot(x,trace2)


    end
end






%%