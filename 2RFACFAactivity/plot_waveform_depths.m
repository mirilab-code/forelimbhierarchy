function plot_waveform_depths(N)
% input the neurons data structure and it will output the waveforms as well as a histogram of the depths

wf = {N.waveform};
d = [N.depth];
d = d-min(d);
d = d/max(d);
d = num2cell(d);


figure;
tiledlayout(1,2);
nexttile;
hold on
lh = cellfun(@(wfm,color) plot(wfm, 'Color',[1-color, 0, 0]), wf,d);
hold off

d = [N.depth];
% u = [N.unit];
w = [N.width];
for i=1:length(d)
    line = lh(i);
    row1 = dataTipTextRow('depth',line.XData*0+d(i));
    row2 = dataTipTextRow('unit',line.XData*0+i);
    row3 = dataTipTextRow('width',line.XData*0+w(i));
    line.DataTipTemplate.DataTipRows(end+1) = row1;
    line.DataTipTemplate.DataTipRows(end+1) = row2;
    line.DataTipTemplate.DataTipRows(end+1) = row3;
end

red = flipud(linspace(0,1,length(d))');
cmap = [red red*0 red*0];
colormap(gca,cmap);
cbar = colorbar;
cbar.Ticks= [0 1];
cbar.TickLabels = {'tip of probe','base of probe'};

nexttile;
histogram([N.depth], 'BinWidth',50);
xlabel('depth from tip')
set(gca,'XDir','reverse');
camroll(-90)


end