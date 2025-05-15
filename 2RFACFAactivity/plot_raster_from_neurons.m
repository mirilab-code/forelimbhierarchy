function plot_raster_from_neurons(N)

units = [N.unit];
figure;
hold on
for i=1:length(units)
    u = units(i);
    scatter(N(units==u).train,zeros(length(N(units==u).train),1)+i, '.k')
end


end