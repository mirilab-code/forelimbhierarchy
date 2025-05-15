function C = find_dense_neighborhood(D,k)
% stupidly inefficient function to find the most dense neighborhood of
% points given a distance matrix. Basically a Vietoris-Rips filtration lol
% Solve the clique problem by interatively binarizing the distance matrix
% filtering by the distance values.

if(size(D,1) < k)
    error('requested clique size is larger than size of matrix');
end

w = triu(D);
filt = sort(w(w~=0));

max_clique = 0;
counter = 1;
while (max_clique<k)
    r = filt(counter);
    A = (D<=r) - eye(size(D));
    
    clique_list = maximalCliques(A);
    [s,~] = cellfun(@size,clique_list);
    max_clique = max(s);
    
    
    
    counter = counter+1;
end
% disp(counter);
max_ind = find(s==max(s));

if (length(max_ind)==1)
    disp(['found unique maximal clique of size ' num2str(max_clique)]);
    C = clique_list(max_ind);
    C = C{1};
else
    disp('non unique maxiaml cliques, returned returned clique with smallest mean distance');
    C = clique_list(max_ind);
    
    clique_size = zeros(length(C),1);
    for i=1:length(C)
        q = D(C{i},C{i});
        clique_size(i) = mean(q(q~=0));
    end
    C = C{find(min(clique_size))};
    
end



% disp(max_ind);
% disp(max_clique);


end