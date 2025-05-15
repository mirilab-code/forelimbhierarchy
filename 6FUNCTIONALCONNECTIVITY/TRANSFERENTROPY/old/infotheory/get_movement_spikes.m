function trn = get_movement_spikes(signal,bintrain,pad)
% from a binary spike train, concatenate only the chunks that align with 1
% on the signal and remove all that are 0. Then pad wit zeros between them
% so when we do the TE with delays you won't have a spike from a different
% movement epoch be calculated into the TE of another one.

% example:
% signal =       [0 0 0 0 1 1 1 1 1 1 1 0 0 0 1 1 1 1 0 0 0]
% bintrain = [0 1 0 0 0 0 1 0 0 1 0 0 1 0 0 1 0 1 1 0 0]
% pad = 4
%
% output = [0 0 1 0 0 1 0 0 0 0 0 0 1 0 1]
%                         |-----|
%                          ^these zeros are from the padding

minT = min(length(signal),length(bintrain));
signal = logical(signal(1:minT));
bintrain = bintrain(1:minT);
padding = zeros(1,pad);

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
    if(i ~= num_epochs) % we don't want to pad the last window?
        chunk = [sub_window padding];
    else
        chunk = sub_window;
    end
    trn = [trn chunk];
end








end