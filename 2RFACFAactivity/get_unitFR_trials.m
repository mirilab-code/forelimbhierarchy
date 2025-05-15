function T = get_unitFR_trials(unit,FR, windows)
% This function gets the firing times of _unit_ in all the windows
% INPUT: unit    := index of a unit in events
%        events  := an nUnitsxTime matrix 
%        windows := an Nx2 matrix with the bounds of the windows we want to get the spike times of
% 
% OUTPUT: T := a nTrialsxWindowLength matrix that should be the same dim as
%               spout trials

nWindows = size(windows,1);
lenWindow = windows(1,2)-windows(1,1)+1;
T = zeros(nWindows,lenWindow);

for i=1:nWindows
    this_window = windows(i,:);
    T(i,:) = FR(unit,this_window(1):this_window(2));
end    

end