function [reach_index, up_time, ramp_time] = doubleThreshold(trace,...
    mean, std, upper_thresh)
%doubleThreshold finds the index where the reach is initiated using
%   a double threshold scheme. The maximal peak is first found and the
%   first most recent threshold cross is marked as the the reach point.
%   trace: EMG signal from -500ms to touch.
%   mean: average of EMG from -500ms to -100ms
%   std: standard deviation of EMG from -500ms to -100ms
%     
    
    %Find the array of all indices where the EMG crosses the .33 std
    %threshold
    cross_index = [];
    for ii = 1:length(trace)-1
        if trace(ii)<mean+7*std
            if trace(ii+1)>=mean+7*std
                cross_index = [cross_index,ii+1];
            end
        end
    end
    %If this threshold cross never happens, we arbitraily assign the cross
    %index to the last index in the trace.
    if (isempty(cross_index))
        cross_index = length(trace)+1;
    end
   
    %Find the array of all indices where the EMG cross the upper
    %threshold
    after_trace = trace(450:end);
    higher_index = [];
    for ii = 1:length(after_trace)-1
        if after_trace(ii)<upper_thresh
            if after_trace(ii+1)>=upper_thresh 
                higher_index = [higher_index,ii+451];
            end
        end
    end
    %If this threshold cross never happens, we arbitraily assign the the
    %the higher index to the length of the trace.
    if (isempty(higher_index))
        higher_index = length(trace)+1;
    end
    
    up_time = higher_index(1);

    %If the cross index is equal to the length of the trace, then there was
    %never a threshold cross. We arbitrily assign the reach index to one
    %minus the length of trace (so that the latency between reach and touch
    %is equal to zero.
    if cross_index == length(trace)+1
        reach_index = length(trace)-1;
        ramp_time = -1;
    elseif higher_index == length(trace)+1
        reach_index = length(trace)-1;
        ramp_time = -1;
    else
        %Find the cross index preceding the peak index.
        diff = -1*(cross_index - higher_index(1));
        for ii = 1:length(diff)
            if diff(ii)<0
                diff(ii) = NaN;
            end
        end

        if isnan(min(diff))
            reach_index = length(trace)-1;
            ramp_time = -1;
        else
            reach_index = higher_index(1)-min(diff);
            ramp_time = min(diff);
        end

    end
end

