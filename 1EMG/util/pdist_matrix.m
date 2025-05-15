function D = pdist_matrix(M)

    nReaches = size(M,3);
    D = zeros(nReaches,nReaches);
    
    for i=1:nReaches
        for j=1:nReaches
            D(i,j) = dtw(M(:,:,i),M(:,:,j));
        end
    end

end