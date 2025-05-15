function bt = train_to_binary(train,maxT)
    %MODIFIED BY DIYA VERY SLIGHTLY

    % this changes a list of spike times to a binary vector with 1's at the
    % spikes
    % NOTE: this does round the spike times into 1ms bins
    times = round(train(train<=maxT));
    times(times==0) = []; %lil confused by this
    bt = zeros(1,maxT+1);
    bt(times) = 1;
    bt = logical(reshape(bt,1,[]));
end