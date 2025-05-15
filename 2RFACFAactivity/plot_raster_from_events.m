function plot_raster_from_events(events,msize,window,plot_title)
    sz = 2;
    if (nargin>=2) 
        sz = msize;
    end
    if (nargin>=3)
        inds = (events(:,2)>=window(1) & events(:,2)<window(end));
        events = events(inds,:);
    end
    hold on
%     scatter(events(:,2),events(:,1), sz, 'k', 'filled');
    scatter(events(:,2),events(:,1), sz, '|k');
    hold off
    xlabel('time (ms)');
    % but it might not be ms so check if your events is in samples or seconds or whatever
    ylabel('neuron index');
    if (nargin>=4)
        title(plot_title);
    end
end