function [train, depth] = events_to_train(events)
    u = unique(events(:,1));
    train = {};
    depth = zeros(size(u));

    for i=1:length(u)
        % get all the spike times for this unit
        train{i} = events(events(:,1)==u(i),2);
        depth(i)= unique(events(events(:,1)==u(i),3));
    end
    train = train';

end