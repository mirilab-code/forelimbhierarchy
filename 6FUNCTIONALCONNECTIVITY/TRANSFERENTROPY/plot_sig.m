function onlysig = plot_sig(pvalmap,units)
    
sig = pvalmap < 0.05;

onlysig = pvalmap;
onlysig(~sig) = nan;

colors = spring(50);
colors = [[1 1 1]; colors];


fignum = randsample(100:200,1);
figure(fignum);
imagesc(onlysig);


totalN = sum(units);
line([units(1),units(1)], [0,totalN], 'Color', 'r');
line( [0,totalN], [units(1),units(1)], 'Color', 'r');

colormap(colors);
h = colorbar;
% set( h, 'YDir', 'reverse' );
caxis([0 0.05]);
ax = gca;
% ax.YDir = 'normal';
xlabel('RFA | CFA');
ylabel('CFA | RFA');
title('row is source, column is target')

figure(fignum+1);
histogram(pvalmap(:));

end