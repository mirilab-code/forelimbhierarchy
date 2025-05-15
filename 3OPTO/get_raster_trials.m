function X = get_raster_trials(train,times,window)
% inputs a spike train along with an array of times and a window of how
% much time to extract around those times. 
% outputs a cell of spike trains.

X = {};

for i=1:length(times)
    t = times(i);
    lb = t+window(1);
    ub = t+window(end);
    trn = train(train>=lb & train<=ub);
    X{i} = trn-t;
end
X = X';

end