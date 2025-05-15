function FR = train_to_firingrate(train)
    % diya mod with progress bar
    % get activity over time
    duration = max(cell2mat(train));
    units = 1:size(train,1);
    
    srate = 1000;                              % Hz 
    min_timevec = 0;                           % sec
    max_timevec = ceil(duration/1000);         % sec
    sigma = 0.001;                              % sec
    peak = 0;
    FR = zeros(length(units),max_timevec*1000+1);
    
    pb = progressbar(' getting firing rate... ');

    for i=1:length(units)
        pb.print(i, length(units));
        if (~isempty(train{i}))
%             disp(i)
            spkvec = train{i}/1000;
            [fr,~,~] = smooth_spikes(spkvec,srate,min_timevec,max_timevec,sigma,peak);
            FR(i,:) = fr'; 
        end
    end

end