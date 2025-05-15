function plot_raster_from_trains(T)
    hold on
    for i=1:length(T)
        x = ones(1,length(T{i}))*i;
        scatter(T{i},x, '|r');
    end
    hold off
    axis tight
end