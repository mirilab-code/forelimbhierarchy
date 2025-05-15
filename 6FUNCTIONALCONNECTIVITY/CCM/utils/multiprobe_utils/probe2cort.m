function cort_neurons = probe2cort(probe)

%from probe data (like neurons_probe1_curated) to getting
%neurons 1400m from surface

%this is somewhat a placeholder function, since I'd imagine there needs to
%be a totally different approach using histology to find relevant neurons.
%for now all this does is some mildly covoluted way to find a surface depth
%and then get all neurons between that and 1400 under. 

%% lets get relevant neurons
depth = [probe.depth];

%lets see drop off in activity near surface
near_surface = depth(depth > 2000); 

[binCounts, binEdges] = histcounts(near_surface, 5);

small_bins = binCounts(binCounts<10);

if length(small_bins) < 1
    surface_depth = binEdges(end);
else

    
    binDiffs = diff(binCounts);
    [~, idxMin] = min(binDiffs);
    
    surface_depth = binEdges(idxMin + 1);
end

bottom = surface_depth-1400; %should be 1400


disp(surface_depth)

cort_neurons = probe(depth > bottom & depth < surface_depth);

%this is extremely coarse, doesn't acct for angle of insertion
%cfa = cfa_probe([cfa_probe.depth] > cfa_surface_depth(end) - 1400); 
%rfa = rfa_probe([rfa_probe.depth] > rfa_surface_depth(end) - 1400);