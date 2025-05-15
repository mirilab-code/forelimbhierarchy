function [trn,slices] = cattrain(signal,bintrain)

minT = min(length(signal),length(bintrain));
signal = logical(signal(1:minT));
bintrain = bintrain(1:minT);
slices = [];

% first get all the windows of spikes and then pad them with 0s
props = regionprops(signal, 'Area');
window_starts = find([0 diff(signal)] == 1);
window_lengths = [props.Area];

if(length(window_starts) ~= length(window_lengths))
%     disp([length(window_starts) length(window_lengths)]);
    window_lengths = window_lengths(2:end);
    % because when we get not_movement from ~movement, the first entry on
    % diff(not_movement) is -1 so we just want to skip over that
end

start_lengths = [window_starts' window_lengths'];
num_epochs = size(start_lengths,1);
trn = [];
for i=1:num_epochs
    sub_window = bintrain(start_lengths(i,1):start_lengths(i,1)+start_lengths(i,2)-1);
    trn = [trn sub_window];
    slices = [slices length(trn)];
end


end