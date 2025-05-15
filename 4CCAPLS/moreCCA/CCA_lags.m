animal = 'ss2';
lags = -30:30;
% lags = 0:0;

R = {};
dCFA = [];
dRFA = [];
VC_CFA = {};
VC_RFA = {};
unweightedRs = [];
weightedRs = [];
for i=lags
    disp(i);
    indx = i-min(lags)+1;
    [r,pcadimCFA,pcadimRFA,vcCFA,vcRFA] = get_CCA_lags(animal,i);
    r = reshape(r,1,[]);
    R{indx} = r;
    dCFA = [dCFA; pcadimCFA];
    dRFA = [dRFA; pcadimRFA];
    VC = [vcCFA vcRFA];
    meanVC = reshape(mean(VC,2),1,[]);    
    weightedR = sum(meanVC .* r);
    
    unweightedRs = [unweightedRs; mean(r)];
    weightedRs = [weightedRs; weightedR];
end
% 
% varcap_total_CFA{end+1} = vcCFA
% varcap_total_RFA{end+1} = vcRFA


%% plot. The lag -30:30 is how much CFA has been shifted.
figure(1);
tiledlayout(2,1);
nexttile;

cmap = winter(length(lags));
title(animal);
hold on
for i=1:length(R)
    plot(R{i},'color',cmap(i,:))
end
hold off
xlabel('canonical variable');
ylabel('R value');
colormap(cmap);
cbar = colorbar;
cbar.Ticks = linspace(0, 1, length(lags(1:5:end)));
cbar.TickLabels = num2cell(lags(1:5:end));
cbar.Label.String = 'CFA lag (ms)';

%
nexttile;
% meanR = cellfun(@mean, R);
plot(lags,weightedRs)
xlabel('lag (ms)');
ylabel('mean of R values');


%%
savename = sprintf('results/%s_CCA.mat',animal);
save(savename, 'R','dCFA','dRFA','unweightedRs','weightedRs','lags','animal');



disp('saved!');
