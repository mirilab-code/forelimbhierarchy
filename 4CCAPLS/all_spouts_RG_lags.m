function ALL_lags = all_spouts_RG_lags(M,spout_bounds,window,lags)

ALL_lags = [];
for i=1:length(lags)
    lag = lags(i);
    lag_bounds = cellfun(@(x) x+lag, spout_bounds, 'UniformOutput', false);
    L = all_spouts_RG(M,lag_bounds,window);
    ALL_lags = cat(3,ALL_lags,L);
end


end