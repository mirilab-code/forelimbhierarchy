% plots rasters of trials in chronological order
% ie, the first trial is on top. RIGHT? 

function plot_raster_from_binary(B) % B is a matrix of [trials x spike train binary]
    hold on
    for t=size(B,1):-1:1 % iterate over trials
        spikes = find(B(t,:));
        x = ones(1,length(spikes))*t;
        scatter(spikes,x, '|k');
    end
    hold off
    axis tight
end