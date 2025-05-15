function FR = train_to_firingrate(train)
    % get activity over time
    duration = max(cell2mat(train)); % in ms
    units = 1:size(train,1);
    
    srate = 1000;                              % Hz 
    min_timevec = 0;                           % sec
    max_timevec = ceil(duration/1000);         % sec
    sigma = 0.01;                              % sec (s.d. of Gaussian)
    peak = 0;
    FR = zeros(length(units),max_timevec*1000+1);

    for i=1:length(units) 
        disp(['Unit ' num2str(i)]) 
        if (~isempty(train{i})) 
            timestamps = train{i}/1000; % convert from spike time in samples (ms) to seconds
            [fr,~,~] = smooth_spikes(timestamps,srate,min_timevec,max_timevec,sigma,peak); 
%             [spkvec,timevec,updatedpeak] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
            FR(i,:) = fr'; 
        end
    end

end