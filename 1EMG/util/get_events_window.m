function E = get_events_window(events,window)
    inds = (events(:,2)>=window(1) & events(:,2)<window(end));
    E = events(inds,:);
end
