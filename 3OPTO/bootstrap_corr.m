function shuffled_corrs = bootstrap_corr(emg,neur,n_iter,df)

n = size(neur,1);
dur = size(neur,2);
minshift = 1000*10; %shift by at least ten seconds

shuffled_corrs = nan(size(emg,1),n, n_iter);
for k=1:n_iter
    fprintf('iter %d out of %d for %s \n',k,n_iter,df);
    neur_shuf = neur;
    r = randi([minshift dur-minshift]);
    neur_shuf = circshift(neur_shuf,r,2);
    
    these_corrs = nan(size(emg,1),n);
    for j=1:size(emg,1)
        these_corrs(j,:) = corr(emg(j,:)',neur_shuf');  
    end

    shuffled_corrs(:,:,k) = these_corrs;
end


end