function ALL = all_spouts_RG(M,spout_bounds,window)

ALL = [];
for i=1:length(spout_bounds)
    sp = spout_bounds{i};
    [reach,grasp] = one_spout_RG(M,sp,window);

    ALL = [ALL reach grasp];

end



end