function T = get_unit_trials(unit,events, windows)
% This function gets the firing times of _unit_ in all the windows
% INPUT: unit    := index of a unit in events
%        events  := an Nx3 matrix with the first column being unit index and the second column being spike time
%        windows := an Nx2 matrix with the bounds of the windows we want to get the spike times of
% 
% OUTPUT: T := a cell with an array of spike times for each window 

nWindows = size(windows,1);
T = cell(nWindows,1);

for i=1:nWindows
    % do the spikes first
    spikes_in_window = (events(:,2)>=windows(i,1) & events(:,2)<windows(i,2));
    E = events(spikes_in_window,:);
    unit_spikes_in_window = E(E(:,1)==unit,:);
    unit_spikes_in_window(:,2) = unit_spikes_in_window(:,2)-windows(i,1);
    
    t = unit_spikes_in_window(:,2);
    
    T{i} = t;
end    

end