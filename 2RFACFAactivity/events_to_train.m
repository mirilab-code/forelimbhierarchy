% EVENTS CANNOT BE A CELL, MUST BE AN (at least) [N x 2] matrix 

function train = events_to_train(events)
    u = unique(events(:,1));
    train = cell(length(u),1);

    for i=1:length(u)
        % get all the spike times for this unit
        unit_spike_times = find(events(:,1)==u(i));
        train{u(i)} = events(unit_spike_times,2);
    end    
end