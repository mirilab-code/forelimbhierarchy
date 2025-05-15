function train = cort_events_to_train(events)
    %cort = events(events(:,3)<=1500, :);
    u = unique(events(:,1));
    train = {};

    for i=1:length(u)
        % get all the spike times for this unit
        unit_spike_times = find(events(:,1)==u(i));
        train{u(i)} = events(unit_spike_times,2);
    end
    train = train';

end