%% plot PLSSVD traces

%% plot traces PLS
lag0 = find(lags==0);
k = 1;
RFAfull = D(k).RFA_RG;
CFAfull = D(k).CFA_RG(:,:,lag0);

dim = 25;

[plsRFA,plsCFA,covs] = pls_svd(RFAfull',CFAfull');

varcapRFA = var(plsRFA);
varcapCFA = var(plsCFA);

varcapRFA = varcapRFA / sum(varcapRFA);
varcapCFA = varcapCFA / sum(varcapCFA);


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
        trace1 = plsRFA(window,comp);
        trace2 = plsCFA(window,comp);

        plot(x,trace1)
        if(j==1)
            covar = covs(comp);
            vc_rfa = varcapRFA(comp);
            vc_cfa = varcapCFA(comp);
            ylabel(sprintf('PLS component %d, \n  covariance: %0.3f \n vcRFA: %0.3f vcCFA: %0.3f ',comp,covar,vc_rfa,vc_cfa))
        end
        yyaxis right
        plot(x,trace2)


    end
end






%%