function FR = singletrain_to_firingrate(train, duration)
    % diya mod with progress bar
    % get activity over time
    %duration = max(cell2mat(train));
    %units = 1:size(train,1);
    
    srate = 1000;                              % Hz 
    min_timevec = 0;                           % sec
    max_timevec = ceil(duration/1000);         % sec
    sigma = 0.001;                              % sec
    peak = 0;
    FR = zeros(1, max_timevec*1000+1);
    
    if (~isempty(train))
%             disp(i)
        spkvec = train/1000;
        [fr,~,~] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
        FR = fr'; 
    end


end