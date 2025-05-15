function P = get_pvals(observed,null,allowed)

P = observed * 0;

for i=1:size(null,1)
    for j=1:size(null,2)
        if(i==j)
            P(i,j) = nan;
        else
            include = squeeze(allowed(i,j,:));
            nulldist = squeeze(null(i,j,include));
            obsval = observed(i,j);

            x = sum(nulldist >= obsval);
            p = (1 + x) / (1 + length(nulldist));
            P(i,j) = p;
        end
    end
end

end