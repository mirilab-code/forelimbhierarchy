
[xc,lags] = spike_train_xcorr(a,b,20);
histogram(xc);

function [xc,lags] = spike_train_xcorr(a,b,maxlag)
    min_diff = min(min(diff(a)),min(diff(b)));
    duration = max(max(a),max(b));
    % now we want to divide [duration] by [min_diff] to get nBins
    nBins = ceil(duration/min_diff);
    
    bin_a = zeros(nBins,1);
    bin_b = zeros(nBins,1);
    
%     bin_a(round(round(a,2)*100)) = 1;
%     bin_b(round(round(b,2)*100)) = 1;
    bin_a(round(a*100)) = 1;
    bin_b(round(b*100)) = 1;

    [xc,lags] = xcorr(bin_a,bin_b,maxlag*100);
    
end