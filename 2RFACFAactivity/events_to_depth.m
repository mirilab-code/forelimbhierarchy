% modified from events_to_train 

function depths = events_to_depth(events)
    u = unique(events(:,1));
    depths = zeros(length(u),1);

    for i=1:length(u)
        % get all the spike times for this unit
        this_neuron = find(events(:,1)==u(i));
        depth = this_neuron(1);
        depths(i) = events(depth,3);
    end
    %depth = depth';

end