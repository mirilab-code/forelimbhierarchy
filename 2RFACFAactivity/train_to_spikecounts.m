function FR = train_to_spikecounts(train)
    % get activity over time
    duration = max(cell2mat(train));
    units = 1:size(train,1);
    
    srate = 1000;                              % Hz 
    min_timevec = 0;                           % sec
    max_timevec = ceil(duration/1000);         % sec
    sigma = 0.01;                              % sec (sd of Gaussian)
    peak = 0;
    FR = zeros(length(units),max_timevec*1000+1);

    for i=1:length(units)
        disp(['Unit ' num2str(i)])
        if (~isempty(train{i}))
            timestamps = train{i}/1000; % convert from spike time in samples to seconds
%             [fr,~,~] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
            [spkvec,spike_count,timevec,~] = smooth_spikes_NK(timestamps,srate,min_timevec,max_timevec,sigma,peak);
            FR(i,:) = spkvec'; 
        end
    end

end