function E = prune_events_by_depth(events, direction, depth)

% inputs an events matrix (with depth of spike in column 3) and removes all
% spikes returns only spike that are above or below a certain threshold.
% 
% INPUTS:  events, direction = above or below, and depth for threshold 
% OUTPUTS: the same type of events matrix but removed all rows with the third
%          column less than/greater than depth
% note: that the depth (third column of events) is smaller at the tip of
%       the probe and larger at the base. so higher numbers are cortex probably

if isequal(direction, 'above')
   inds = find(events(:,3) >= depth); 
end
if isequal(direction, 'below')
   inds = find(events(:,3) <= depth);
end

E = events(inds,:);

end


