function FR = train_to_FR(train,units)
    % get activity over time
    train = reshape(train,[],1);
    duration = max(cell2mat(train));
    
    srate = 1000;                              % Hz 
    min_timevec = 0;                           % sec
    max_timevec = ceil(duration/1000);         % sec
    sigma = 0.01;                              % sec
    peak = 0;

    FR = zeros(length(units),max_timevec*1000+1);

    for i=1:length(units)
        if (~isempty(train{units(i)}))
            spkvec = train{units(i)}/1000;
            [fr,~,~] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
            FR(i,:) = fr'; 
        end
    end

end