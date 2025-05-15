function P = get_pvals_HFS(observed,null,allowed)

P = zeros(length(null),1);

for i=1:length(null)
  
            include = allowed {i};
            nulldist = null{i};
            obsval = observed(i,3);

            x = sum(nulldist >= obsval);
            p = (1 + x) / (1 + length(nulldist));
            P(i) = p;
        
    
end

end